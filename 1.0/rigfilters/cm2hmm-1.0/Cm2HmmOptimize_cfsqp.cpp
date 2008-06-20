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
#include "stdafx.h"

#include "UseDebugNew.h"

#include "cmzasha.h"

// optimize HMM using cfsqp.c functions

//#ifdef _MSC_VER
#define CFSQP_STDC
//#endif
extern "C" {
#include <cfsqpusr.h>
}

class SolverClass_cfsqp {
protected:

	struct CookieData {
		ObjectiveFunc *objectiveFunc;
		std::vector<const Inequality *> constraints;

		std::vector<double> objFuncCachedProblemVars;
		double objFuncCachedResult;
		std::vector<double> objGradientCachedProblemVars;
		std::vector<double> objGradientCachedResult;
		SolverWrapper::MessageReceiver *messageReceiver;
	};

	static void CopyVars_Cfsqp2Vector (std::vector<double>& vars,const double *x,int numVars);
	static void CopyVars_Vector2Cfsqp (double *x,const std::vector<double>& vars,int numVars);

	static void ObjectiveFunction(int nparam,int j,double *x,double *fj,void *voidCookieData);
	static void ObjectiveFunction_Gradient(int nparam,int j,double *x,double *gradfj,void (*dummy)(int,int,double *,double *,void *),void *voidCookieData);
	static void ConstraintFunction(int nparam,int j,double *x,double *gj,void *voidCookieData);
	static void ConstraintFunction_Gradient (int nparam,int j,double *x,double *gradgj,void (*dummy)(int,int,double *,double *,void *),void *voidCookieData);

	// for B,C, see p. 18 of the user manual, describing the 'mode' input parameter
	int B,C;

public:
	SolverClass_cfsqp (int B_,int C_);
	~SolverClass_cfsqp ();
	std::vector<double> Solve(ObjectiveFunc& objectiveFunc,const std::vector<double>& inputProblemVars,bool importantBoundsAreSet,double importantLowerBoundAllVars,double importantUppderBoundAllVars,SolverWrapper::MessageReceiver *messageReceiver) const;
};
SolverClass_cfsqp::SolverClass_cfsqp (int B_,int C_)
{
	B=B_;
	C=C_;
}
SolverClass_cfsqp::~SolverClass_cfsqp ()
{
}
std::vector<double> SolverClass_cfsqp::Solve(ObjectiveFunc& objectiveFunc,const std::vector<double>& inputProblemVars,bool importantBoundsAreSet,double importantLowerBoundAllVars,double importantUpperBoundAllVars,SolverWrapper::MessageReceiver *messageReceiver) const
{
	std::vector<double> problemVars=inputProblemVars;

	int ndim=objectiveFunc.GetNumProblemVars();

	CookieData cookieData;
	cookieData.objectiveFunc=&objectiveFunc;
	cookieData.messageReceiver=messageReceiver;

	// set up linear inequalities constraints
	const InequalityList& srcInequalityList=objectiveFunc.GetInequalityList();
	cookieData.constraints.reserve(srcInequalityList.size());
	for (InequalityList::const_iterator ineqIter=srcInequalityList.begin(); ineqIter!=srcInequalityList.end(); ineqIter++) {
		cookieData.constraints.push_back(&(*ineqIter));
	}

	// code copied from sampl1.c from the CFSQP distribution

	int nparam,nf,nineq,neq,mode,iprint,miter,neqn,nineqn,
		ncsrl,ncsrn,nfsr,mesh_pts[1],inform;
	double bigbnd,eps,epsneq,udelta;
	double *x,*bl,*bu,*f,*g,*lambda;

	int A=0; // we'll always want this
	mode=C*100 + B*10 + A;
	iprint=0;
	miter=500;  
	bigbnd=1.e10;
	eps=1.e-7;
	epsneq=0.e0;
	udelta=0.e0;
	nparam=ndim;
	nf=1;
	neqn=0;
	nineqn=0;
	nineq=(int)(srcInequalityList.size());
	neq=0;
	ncsrl=ncsrn=nfsr=mesh_pts[0]=0;
	bl=(double *)calloc(nparam,sizeof(double));
	bu=(double *)calloc(nparam,sizeof(double));
	x=(double *)calloc(nparam,sizeof(double));
	f=(double *)calloc(nf,sizeof(double));
	g=(double *)calloc(nineq+neq,sizeof(double));
	lambda=(double *)calloc(nineq+neq+nf+nparam,sizeof(double));

	if (importantBoundsAreSet) {
		// implement important bounds
		for (int i=0; i<ndim; i++) {
			bl[i]=importantLowerBoundAllVars;
			bu[i]=importantUpperBoundAllVars;
		}
	}
	else {
		// sets upper & lower bounds on variables to -/+ infinity (i.e., they're unbounded)
		for (int i=0; i<ndim; i++) {
			bl[i]=-bigbnd;
			bu[i]=+bigbnd;
		}
	}

	CopyVars_Vector2Cfsqp(x,problemVars,objectiveFunc.GetNumProblemVars());

	cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
		mode,iprint,miter,&inform,bigbnd,eps,epsneq,udelta,bl,bu,x,
		f,g,lambda,ObjectiveFunction,ConstraintFunction,ObjectiveFunction_Gradient,ConstraintFunction_Gradient,&cookieData);
	//printf("(cfsqp.c) inform==%d\n",inform);
	if (inform!=0) {
		printf("WARNING: (cfsqp.c) inform==%d\n",inform);
		/*
		throw SimpleStringException("WARNING: (cfsqp.c) inform==%d\n",inform);
		*/
	}

	// retrieve solution
	CopyVars_Cfsqp2Vector(problemVars,x,objectiveFunc.GetNumProblemVars());

	// free arrays used for cfsqp
	free(bl);
	free(bu);
	free(x);
	free(f);
	free(g);
	free(lambda);

	return problemVars;
}
void SolverClass_cfsqp::ObjectiveFunction(int nparam,int j,double *x,double *fj,void *voidCookieData)
{
	CookieData *cookieData=(CookieData *)voidCookieData;
	assert(cookieData->objectiveFunc->GetNumProblemVars()==nparam);
	assert(j==1); // only one objective function, and j is 1-based

	std::vector<double> problemVars;
	CopyVars_Cfsqp2Vector(problemVars,x,cookieData->objectiveFunc->GetNumProblemVars());

	if (problemVars==cookieData->objFuncCachedProblemVars) {
		*fj=cookieData->objFuncCachedResult;
	}
	else {
		double fx;
		std::vector<double> gradient;
		vector2d<double> hessian;
		cookieData->objectiveFunc->Eval (fx,gradient,hessian,problemVars,false,false);
		*fj=fx;

		if (cookieData->messageReceiver!=NULL) {
			cookieData->messageReceiver->EvaluatedObjectiveFunc (fx,problemVars);
		}

		cookieData->objFuncCachedProblemVars=problemVars;
		cookieData->objFuncCachedResult=fx;
#if 0
		unsigned int *pui=(unsigned int *)&fx;
		printf("fx=%lg (%lx%lx)\n",fx,pui[1],pui[0]);
#endif
	}
}
void SolverClass_cfsqp::ObjectiveFunction_Gradient(int nparam,int j,double *x,double *gradfj,void (*dummy)(int,int,double *,double *,void *),void *voidCookieData)
{
	CookieData *cookieData=(CookieData *)voidCookieData;
	assert(cookieData->objectiveFunc->GetNumProblemVars()==nparam);
	assert(j==1); // only one objective function, and j is 1-based

	std::vector<double> problemVars;
	CopyVars_Cfsqp2Vector(problemVars,x,cookieData->objectiveFunc->GetNumProblemVars());

	double fx;
	std::vector<double> gradient;
	if (problemVars==cookieData->objGradientCachedProblemVars) {
		gradient=cookieData->objGradientCachedResult;
	}
	else {
		vector2d<double> hessian;
		cookieData->objectiveFunc->Eval (fx,gradient,hessian,problemVars,false);

		cookieData->objGradientCachedProblemVars=problemVars;
		cookieData->objGradientCachedResult=gradient;
	}

	for (int i=0; i<nparam; i++) {
		gradfj[i]=gradient[i];
	}

#if 0
	printf("g=");
	for (int i=0; i<nparam; i++) {
		printf(" %lg",gradient[i]);
	}
	printf("\n");
#endif
}
void SolverClass_cfsqp::ConstraintFunction(int nparam,int j,double *x,double *gj,void *voidCookieData)
{
	j--; // csqfp has this 1-based, for some reason
	CookieData *cookieData=(CookieData *)voidCookieData;
	assert(cookieData->objectiveFunc->GetNumProblemVars()==nparam);
	const Inequality& ineq=*(cookieData->constraints[j]);

	std::vector<double> problemVars;
	CopyVars_Cfsqp2Vector(problemVars,x,cookieData->objectiveFunc->GetNumProblemVars());

	// we have an inequality of the form
	// x_i+x_j+x_k >= 5.
	// to make cfsqp happy, it looks like we'd like to change it to a less than, with 0 on the rhs, i.e.
	// 5-x_i-x_j-x_k <= 0
	double value=ineq.rhs;
	if (ineq.inequalityType==IneqType_GE) {
		value -= 3e-6; // numerical issues
	}
	for (std::list<InequalityTerm>::const_iterator termIter=ineq.lhs.begin(); termIter!=ineq.lhs.end(); termIter++) {
		value -= problemVars[termIter->variableNum];
	}
	if (ineq.inequalityType==IneqType_Less) {
		value = -value; // reverse sign
	}

	*gj=value;
}
void SolverClass_cfsqp::ConstraintFunction_Gradient (int nparam,int j,double *x,double *gradgj,void (*dummy)(int,int,double *,double *,void *),void *voidCookieData)
{
	j--; // csqfp has this 1-based, for some reason
	CookieData *cookieData=(CookieData *)voidCookieData;
	assert(cookieData->objectiveFunc->GetNumProblemVars()==nparam);
	const Inequality& ineq=*(cookieData->constraints[j]);

	std::vector<double> problemVars;
	CopyVars_Cfsqp2Vector(problemVars,x,cookieData->objectiveFunc->GetNumProblemVars());

	// clear gradient
	for (int i=0; i<nparam; i++) {
		gradgj[i]=0;
	}
	for (std::list<InequalityTerm>::const_iterator termIter=ineq.lhs.begin(); termIter!=ineq.lhs.end(); termIter++) {
		switch (ineq.inequalityType) {
			case IneqType_GE:
				gradgj[termIter->variableNum]=-1;
				break;
			case IneqType_Less:
				gradgj[termIter->variableNum]=+1;
				break;
			default:
				assert(false);
		}
	}
}
void SolverClass_cfsqp::CopyVars_Cfsqp2Vector (std::vector<double>& vars,const double *x,int numVars)
{
	vars.resize(numVars);
	for (int i=0; i<numVars; i++) {
		vars[i]=x[i];
	}
}
void SolverClass_cfsqp::CopyVars_Vector2Cfsqp (double *x,const std::vector<double>& vars,int numVars)
{
	assert(vars.size()==(size_t)numVars);
	for (int i=0; i<numVars; i++) {
		x[i]=vars[i];
	}
}

class SolverWrapper_cfsqp : public SolverWrapper {
protected:
	int B,C;
public:
	SolverWrapper_cfsqp (int B_,int C_);
	~SolverWrapper_cfsqp ();

	std::vector<double> /* optimal problem vars */ Solve (ObjectiveFunc *objectiveFunc,const std::vector<double>& inputProblemVars,double maxVariableMagnitudeForUpperLowerBounds,bool importantBoundsAreSet,double importantLowerBoundAllVars,double importantUppderBoundAllVars,MessageReceiver *messageReceiver);
};
SolverWrapper_cfsqp::SolverWrapper_cfsqp (int B_,int C_)
{
	B=B_;
	C=C_;
}
SolverWrapper_cfsqp::~SolverWrapper_cfsqp ()
{
}
std::vector<double> SolverWrapper_cfsqp::Solve (ObjectiveFunc *objectiveFunc,const std::vector<double>& inputProblemVars,double maxVariableMagnitudeForUpperLowerBounds,bool importantBoundsAreSet,double importantLowerBoundAllVars,double importantUppderBoundAllVars,MessageReceiver *messageReceiver)
{
	SolverClass_cfsqp solver(B,C);
	return solver.Solve(*objectiveFunc,inputProblemVars,importantBoundsAreSet,importantLowerBoundAllVars,importantUppderBoundAllVars,messageReceiver);
}

SolverWrapper *NewSolverWrapper_cfsqp (int B,int C)
{
	// for B,C, see p. 18 of the user manual, describing the 'mode' input parameter
	assert(B==0 || B==1);
	assert(C==1 || C==2);
	return new SolverWrapper_cfsqp(B,C);
}
