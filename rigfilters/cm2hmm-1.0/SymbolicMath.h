/*
This file copyright (c) 2004, Zasha Weinberg
All rights reserved.

Redistribution and use in source and binary forms, 
with or without modification, are permitted 
provided that the following conditions are met:

- Redistributions of source code must retain the 
above copyright notice, this list of conditions 
and the following disclaimer. 
- Redistributions in binary form must reproduce 
the above copyright notice, this list of 
conditions and the following disclaimer in the 
documentation and/or other materials provided 
with the distribution. 
- Neither the name of the University of Washington 
nor the names of its contributors may be used to 
endorse or promote products derived from this 
software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// semi-generic class to do symbolic math, such as is used by the symbolic differentiation of the infinite-length forward alg
class SymbolicMath {
	friend class Expression;

public:
	class ExpressionNode {
		// because expanding the expression is exponential, we must cache values, whether these values are computed by Eval, Derivative or DoubleDerivative
	private: // only ExpressionNode will worry about this

		bool isValueClear;
		bool isVisited; // used by various walking functions
		struct DerivativeValue {
			int wrtVarNum;
			double value;
			bool isValid;
		};
		bool isEvalValueValid;
		double evalValue;
		DerivativeValue derivativeValue[2]; // needs 2, for calculation of DoubleDerivative
		bool isDoubleDerivativeValueValid;
		double doubleDerivativeValue;

		// I've got to reference counting, since (1) in building the expression up, many subtrees get discarded, i.e. there is no path to them from the expression root, so we can't just delete them from the root, (2) C++ doesn't provide garbage collection, (3) having objects add themselves to a list isn't threadsafe & I don't want to be seriously un-threadsafe, (4) adding objects to a per-SymbolicProbVariableMath list requires giving each Expression a pointer to SymbolicProbVariableMath, which is inconvenient in the ScanHmm code, which assumes it's dealing with something like a number.
		int refCount;
	protected:
		bool IsValueClear (void) const;
		virtual double ActualEval (const vector<double>& globalVars) = 0;
		virtual double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo) = 0;
		virtual double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2) = 0; // evaluate derivative wrt to var1, then var2
		void Internal_DumpSubtreeEvalCCode (FILE *out);
	public:
		ExpressionNode (void);
		virtual ~ExpressionNode ();
		// NOTE: before calling Eval, Derivative or DoubleDerivative, you must call ClearValue on the root
		double Eval (const vector<double>& globalVars);
		double Derivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double DoubleDerivative (const vector<double>& globalVars,int var1,int var2);

		virtual bool IsConst (void); // including children.  default: if no children, then false, else explore children and return the AND of the children
		virtual int GetNumChildren (void); // default: no children
		virtual ExpressionNode *GetChild (int child);
		virtual double ToConstDouble (void); // throws if not overriden; throws if !IsConst
		virtual void ClearValue (void); // default implementation assumes no children

		void IncRef (void);
		void DecRef (void); // calls delete this, if necessary

		// WARNING: if the expression is even moderately large, this is very bad because (1) it's exponential for a DAG, and (2) everything's fit onto 1 line.  It's only really useful for debugging.
		virtual void DumpExpandedOneLine (FILE *out); // default is to throw an exception
		// prints one line of C code to eval this; assumes nodes correspond to variables of the form t# where # is the object's this ptr in hex.  The caller is responsible for ensuring that the calls happen in a feasible order.  The code should be compatible with C.  the globalVars are in an array of doubles.  varToDifferentiateTo and var1,var2 are defined as before
		virtual void DumpEvalCCode (FILE *out); // default: throws 'not implemented'
		// dumps the subtree rooted here into CCode
		void DumpSubtreeEvalCCode (FILE *out);

		// poor man's run-time type identification (I don't want to have to enable RTTI for everything, just so this class will work, so I'll use the poor man's approach)
		virtual bool Is_SumOfConstantTimesExpression (void) const;
		virtual bool Is_BinaryMult (void) const;
		virtual bool Is_LiteralConst (void) const;
	};
	typedef std::list<ExpressionNode *> ExpressionNodeList;

	class ExpressionNode_Null : public ExpressionNode { // dummy node, for default constructor
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Null (void);
		~ExpressionNode_Null ();
	};
	class ExpressionNode_Const : public ExpressionNode {
	protected:
		double x;
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Const (double t);
		~ExpressionNode_Const ();
		bool IsConst (void);
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
		bool Is_LiteralConst (void) const;
	};
	class ExpressionNode_VarPow2 : public ExpressionNode {
	protected:
		int varNum;
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_VarPow2 (int varNum_);
		~ExpressionNode_VarPow2 ();
		void DumpExpandedOneLine (FILE *out);
	};
	class ExpressionNode_Var : public ExpressionNode {
	protected:
		int varNum;
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Var (int varNum_);
		~ExpressionNode_Var ();
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_BinaryOp : public ExpressionNode {
	protected:
		ExpressionNode *f,*g;
	public:
		ExpressionNode_BinaryOp (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_BinaryOp ();
		void ClearValue (void);
		int GetNumChildren (void);
		ExpressionNode *GetChild (int child);
	};
	class ExpressionNode_Add : public ExpressionNode_BinaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Add (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_Add ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_Minus : public ExpressionNode_BinaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Minus (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_Minus ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_Mult : public ExpressionNode_BinaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Mult (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_Mult ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
		bool Is_BinaryMult (void) const;
	};
	class ExpressionNode_Div : public ExpressionNode_BinaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Div (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_Div ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_Pow : public ExpressionNode_BinaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Pow (ExpressionNode *f_,ExpressionNode *g_);
		~ExpressionNode_Pow ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_UnaryOp : public ExpressionNode {
	protected:
		ExpressionNode *f;
	public:
		ExpressionNode_UnaryOp (ExpressionNode *f_);
		~ExpressionNode_UnaryOp ();
		void ClearValue (void);
		int GetNumChildren (void);
		ExpressionNode *GetChild (int child);
	};
	class ExpressionNode_Log2 : public ExpressionNode_UnaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Log2 (ExpressionNode *f_);
		~ExpressionNode_Log2 ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_Sqrt : public ExpressionNode_UnaryOp {
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Sqrt (ExpressionNode *f_);
		~ExpressionNode_Sqrt ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	class ExpressionNode_MultiParamOp : public ExpressionNode {
	protected:
		typedef vector<ExpressionNode *> ExpressionNodeList; // for quicker GetChild calls
		ExpressionNodeList expressionNodeList;
	public:
		ExpressionNode_MultiParamOp (void);
		~ExpressionNode_MultiParamOp ();

		void AppendParam (ExpressionNode *f);
		void ClearValue (void);
		int GetNumChildren (void);
		ExpressionNode *GetChild (int child);
	};
	class ExpressionNode_Summation : public ExpressionNode_MultiParamOp { // could do this with ExpressionNode_Add, but with many terms to add, it can overflow the stack, and trying to balance the tree is a major hassle
	protected:
		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
	public:
		ExpressionNode_Summation (void);
		~ExpressionNode_Summation ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
	};
	// this class is used to add common expressions with different constants, as we get in the inf-len forward alg
	class ExpressionNode_SumOfConstantTimesExpression : public ExpressionNode {
	protected:
		struct Term {
			double factor;
			ExpressionNode *expressionNode;
			inline bool operator < (const Term& t) const {
				return expressionNode < t.expressionNode; // sort by sub-expression identifiers (i.e. pointers)
			}
		};
		typedef vector<Term> TermList;
		TermList termList;

		double ActualEval (const vector<double>& globalVars);
		double ActualDerivative (const vector<double>& globalVars,int varToDifferentiateTo);
		double ActualDoubleDerivative (const vector<double>& globalVars,int var1,int var2);
		void ExtractTerms(ExpressionNode *f,double factorSoFar=1);
		void CombineLikeTerms(void);
	public:
		// externally behaves like Plus
		ExpressionNode_SumOfConstantTimesExpression (ExpressionNode *f,ExpressionNode *g);
		~ExpressionNode_SumOfConstantTimesExpression ();
		double ToConstDouble (void);
		void DumpExpandedOneLine (FILE *out);
		void DumpEvalCCode (FILE *out);
		void ClearValue (void);
		int GetNumChildren (void);
		ExpressionNode *GetChild (int child);
		bool Is_SumOfConstantTimesExpression (void) const;
	};
private: // I'd like to protect derived classes from this -- they should call SetExpressionNode
	ExpressionNode *rootExpressionNode;
	void DeleteAllExpressionNode(void);
protected:
	// deferred: void SetRootExpression (const Expression& rootExpression);
	void SetRootExpression (ExpressionNode *rootExpressionNode_);
	double Eval (const vector<double>& problemVars);
	double Derivative (const vector<double>& problemVars,int problemVarToDifferentiateTo);
	double DoubleDerivative (const vector<double>& problemVars,int problemVarToDifferentiateTo_1,int problemVarToDifferentiateTo_2);

public:
	SymbolicMath (void);
	~SymbolicMath ();

	void Eval (int numVars,double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient=true);

	// the actual class that has operations performed on it
	class Expression {
	protected:

		/// WARNING: enableGroupingCommonTerms should probably be set to false for most problems; it only works with something like the inf-len fwd alg, where there are many like terms that can be added together.  Without this, it might take O(n^2) for 'n' terms with a tree containing lots of adds

		static bool enableGroupingCommonTerms; // for now, default to true


		static bool enableBinaryOpCache; // just uses lots of RAM & doesn't help that much
		ExpressionNode *expressionNode;
		static ExpressionNode_Null nullExpressionNode;
		bool HasSymmetricAnnihilator(const Expression& t,double annihilator);
		bool HasSymmetricIdentityConst(const Expression& t,double identity);
		void CheckForConst(void);
		static ExpressionNode *CreateConst (double x);

		class AutoDecExpressionNode {
		protected:
			ExpressionNode *p;
			void DecRef (void) {
				if (p!=NULL) {
					p->DecRef();
				}
			}
		public:
			AutoDecExpressionNode (void) { p=NULL; }
			void operator = (const AutoDecExpressionNode& t) { DecRef(); p=t.p; p->IncRef(); }
			AutoDecExpressionNode (const AutoDecExpressionNode& t) { p=NULL; *this=t; }
			void operator = (ExpressionNode *t) { DecRef(); p=t; p->IncRef(); }
			AutoDecExpressionNode (ExpressionNode *t) { p=NULL; *this=t; }
			operator ExpressionNode * () { return p; }
			~AutoDecExpressionNode () {
				DecRef();
			}
		};

		typedef std::map<double,AutoDecExpressionNode> ConstMap;
		static ConstMap constMap;

		enum OpType {
			OpType_Mult,OpType_Add
		};
		struct BinaryOpDef {
			int opType;
			ExpressionNode *f,*g;
			bool operator < (const BinaryOpDef& t) const {
				if (opType!=t.opType) {
					return opType<t.opType;
				}
				if (f!=t.f) {
					return f<t.f;
				}
				return g<t.g;
			}
		};
		void PostprocessSymmetricBinaryOpForCache(BinaryOpDef thisOp);
		typedef std::map<BinaryOpDef,AutoDecExpressionNode> BinaryOpMap;
		static BinaryOpMap binaryOpMap;
	public:
		void operator = (const Expression& t);
		void operator = (double t);
		Expression (void);
		Expression (const Expression& t);
		Expression (double t);
		Expression (ExpressionNode *t);
		void operator *= (const Expression& t);
		void operator += (const Expression& t);
		void operator -= (const Expression& t);
		void operator /= (const Expression& t);
		static Expression ExpressionOfVarPow2 (int var); // sets to 2^(var)
		static Expression ExpressionOfVar (int var); // sets to var
		static Expression Log2 (const Expression& t);
		static Expression Sqrt (const Expression& t);
		static Expression Pow (const Expression& mantissa,int exponent);
		static Expression Pow (const Expression& mantissa,const Expression& exponent);
		ExpressionNode *GetExpressionNode (void);
		~Expression ();

		double Eval (const vector<double>& problemVars);
		void DumpExpandedOneLine (FILE *out); // see warning under same-named func in ExpressionNode
		void DumpEvalCCode (FILE *out);

		// can be good to periodically clear the cache, to free it of useless expressions.  I've now disabled the general common sub-expr cache, so that's irrelevant, but clearing consts is good, since many random constants don't get re-used.  (Actually, I should probably move to a thingy where a const object frees itself from the cache once it's done with, but I'm too lazy.)
		static void ClearCommonSubExpressionsCache (void);
		static void ClearConstCache (void);

		// throws exception if it's not of type const
		double ToConstDouble (void);
	};
protected:
	void SetRootExpression (Expression& rootExpression);
public:
	SymbolicMath (Expression expression);

	class MultiParamExpression {
	protected:
		ExpressionNode_MultiParamOp *expressionNode;
		MultiParamExpression (ExpressionNode_MultiParamOp *t);
	public:
		~MultiParamExpression ();
		Expression ToExpression ();
		void AppendParam (Expression& t);
	};
	class SummationExpression : public MultiParamExpression {
	public:
		SummationExpression ();
	};
};
SymbolicMath::Expression operator + (const SymbolicMath::Expression& x,const SymbolicMath::Expression& y);
SymbolicMath::Expression operator - (const SymbolicMath::Expression& x,const SymbolicMath::Expression& y);
SymbolicMath::Expression operator * (const SymbolicMath::Expression& x,const SymbolicMath::Expression& y);
SymbolicMath::Expression operator / (const SymbolicMath::Expression& x,const SymbolicMath::Expression& y);
