//---------------------------------------------------------------------------
#pragma hdrstop

#include "TwoSampleZhangHomodenTests.h"
#pragma package(smart_init)
#include "define.h"
#include "Estimators/Estimator.h"
#include "MathLib/IsMathLib.h"
#include "ParametersLib/ParametersLib.h"
#include "CriteriaLib/criterion_lib.h"
#include "ManagerLib/TestManager.h"

#include <algorithm>

// Структура для вычисления рангов

	struct elem
	{
	 public:
	 slongint ind;
	 slongint nSapmle;
	 slongfloat y;
	 elem::elem(){ind=0;y=0.;nSapmle=0;}
	};

	bool ssort (const elem &e1,const elem &e2)
	{
	  return (e1.y < e2.y);
	}
//---------------------------------------------------------------------------
/// Критерий "Zk"
 Criterion_ZhangHomodenZk::Criterion_ZhangHomodenZk()
   {
		_info.SetName("Критерий Жанга Zk (2выборки)");
		_info.SetFileName("Критерий Жанга Zk (2выборки)");
		_info.SetNote("Критерий Жанга Zk (2выборки)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_ZhangHomodenZk::~Criterion_ZhangHomodenZk(){}

   void Criterion_ZhangHomodenZk::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_ZhangHomodenZk::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());
	   sint m=(*SampVec)[0]->n();  //x
	   sint n=(*SampVec)[1]->n();  //y
	   sint size = m+n;

	// Количество выборок
	sint samples_number = (*SampVec).size();
	vector <sint> j_i(samples_number, 0);

	slongfloat  Fk, Fik;
	slongfloat Zk = -1e100,sum=0.0;



	vector<vector<sint> > R(samples_number);
	vector <elem> x(0);

	//вычислим ранги Rij
	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
			for (sint j = 0; j < ni; j++)
		{
			elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
			#pragma omp critical
		   {
			x.push_back(el);
		   }

		}
	}

	sort( x.begin(), x.end(), ssort);
	// получаем ранг каждого наблюдения
	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
		vector < sfloat> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
			R[i][j]=buf[j] ;

	}

	for(sint k = 1; k <= size; k++)
	{
		sum=0.0;
		for(sint i = 0; i < samples_number; i++)
		{
		if(k==R[i][j_i[i]])
			{
				j_i[i]++;
				Fik = (j_i[i]-0.5)/sfloat((*SampVec)[i]->n());
			}
			else
				Fik = j_i[i]/sfloat((*SampVec)[i]->n());
			Fk = (k-0.5)/sfloat(size);
			if(  fabs(Fik)<1e-12)
				sum += ((*SampVec)[i]->n())*( (1-Fik)*log((1-Fik)/(1-Fk)));
			 else
				if(  fabs(1-Fik)<1e-12)
						sum += ((*SampVec)[i]->n())*(Fik * log(Fik/Fk)); 
				else
					sum += ((*SampVec)[i]->n())*(Fik * log(Fik/Fk) + (1-Fik)*log((1-Fik)/(1-Fk)));
		}

		Zk = dmax(Zk, sum);
	}
	   return Zk;
   }

   CriterionBase* Criterion_ZhangHomodenZk::Copy()
   {
	  return (new Criterion_ZhangHomodenZk());
   }

   Distribution * Criterion_ZhangHomodenZk::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
  }

 //----------------------------------------------------------------------------------------
 //---------------------------------------------------------------------------
/// Критерий "Za"
  Criterion_ZhangHomodenZa::Criterion_ZhangHomodenZa()
   {
		_info.SetName("Критерий Жанга Za (2выборки)");
		_info.SetFileName("Критерий Жанга Za (2выборки)");
		_info.SetNote("Критерий Жанга Za (2выборки)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_ZhangHomodenZa::~Criterion_ZhangHomodenZa(){}

   void Criterion_ZhangHomodenZa::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_ZhangHomodenZa::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());
	   sint m=(*SampVec)[0]->n();  //x
	   sint n=(*SampVec)[1]->n();  //y
	   sint size = m+n;

	sint samples_number = (*SampVec).size();
	vector <slongfloat> j_i(samples_number, 0.0);

	slongfloat  Fik;
	slongfloat Za = 0.0,sum=0.0;

	vector<vector<sint> > R(samples_number);
	vector <elem> x(0);
	//вычислим ранги Rij

	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
			for (sint j = 0; j < ni; j++)
		{
			elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
			#pragma omp critical
		   {
			x.push_back(el);
		   }
		}
	}
	sort( x.begin(), x.end(), ssort);

// получаем ранг каждого наблюдения
	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
		vector < sint> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
					R[i][j]=buf[j] ;
		buf.clear();
	}

	for(sint k = 1; k <= size; k++)
		for(sint i = 0; i < samples_number; i++)
		{
		  slongfloat ni= (*SampVec)[i]->n();
			if(k==R[i][j_i[i]])
			{
				j_i[i]=j_i[i]+1;
				Fik = (j_i[i]-0.5)/ni;
			}
			else
				Fik = j_i[i]/ni;


			if(  fabs(Fik)<1e-12 || Fik==0)
				Za += ni*(1.-Fik)*log(1.0-Fik)/((k-0.5)*(size-k+0.5));
			 else
				if(  fabs(1-Fik)<1e-12 || Fik==1 )
					Za += ni * Fik * log(Fik)/((k-0.5)*(size-k+0.5));
				else
					if(k==size)
						Za += ni *(Fik * log(Fik) + (1.0-Fik)*log(1.0-Fik))/sfloat((k-0.5)*(size-k+0.5));
					else
						Za += ni *(Fik * log(Fik) + (1.0-Fik)*log(1.0-Fik))/sfloat((k-0.5)*(size-k+0.5));
		}

	return -Za;
   }

   CriterionBase* Criterion_ZhangHomodenZa::Copy()
   {
	  return (new Criterion_ZhangHomodenZa());
   }

   Distribution * Criterion_ZhangHomodenZa::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
   }
   //----------------------------------------------------------------------------------------
//// Критерий "Zc"
 Criterion_ZhangHomodenZc::Criterion_ZhangHomodenZc()
   {
		_info.SetName("Критерий Жанга Zc (2выборки)");
		_info.SetFileName("Критерий Жанга Zc (2выборки)");
		_info.SetNote("Критерий Жанга Zc (2выборки)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_ZhangHomodenZc::~Criterion_ZhangHomodenZc(){}

   void Criterion_ZhangHomodenZc::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_ZhangHomodenZc::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());
	   sint m=(*SampVec)[0]->n();  //x
	   sint n=(*SampVec)[1]->n();  //y
	   sint size = m+n;

	sfloat Zc = 0.0;
	vector<vector<sfloat> > R(2);
	vector <elem> x(0);
	//Вычислим ранги Rij
   // x.reserve(size);
	for(sint i = 0; i < 2; i++)
	{
		sint ni= (*SampVec)[i]->n();
		R[i].reserve(ni);
		for (sint j = 0; j < ni; j++)
		{
		 	elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
			#pragma omp critical
		   {
			x.push_back(el);
		   }
		}
	}
	 /*	x[0].y=2.0;
	x[1].y=4.0;
	x[2].y=18.0;
	x[3].y=9.0;
	x[4].y=10.0;
	x[5].y=25.0;     */
	sort( x.begin(), x.end(), ssort);

	// получаем ранг каждого наблюдения
	for(sint i = 0; i < 2; i++)
	{
		sint ni= (*SampVec)[i]->n();
		vector < sfloat> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
					R[i][j]=buf[j] ;
		buf.clear();
	}
	for(sample_n_t i = 0; i < 2; i++)
	{
		sfloat ni=(*SampVec)[i]->n();
		for(sample_n_t j = 0; j < ni; j++)
			Zc +=  log(ni/(j+1-0.5)-1.)*log(sfloat(size)/(R[i][j]-0.5)-1.);

	}
	return Zc/sfloat(size);
   }

   CriterionBase* Criterion_ZhangHomodenZc::Copy()
   {
	  return (new Criterion_ZhangHomodenZc());
   }

   Distribution * Criterion_ZhangHomodenZc::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
}
   //----------------------------------------------------------------------------------------

//---------------------------------------------------------------------------
/// K-выборочный Критерий "Zk"
 Criterion_KZhangHomodenZk::Criterion_KZhangHomodenZk()
   {
		_info.SetName("Критерий Жанга Zk (К выборок)");
		_info.SetFileName("Критерий Жанга Zk (К выборок)");
		_info.SetNote("Критерий Жанга Zk (К выборок)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_KZhangHomodenZk::~Criterion_KZhangHomodenZk(){}

   void Criterion_KZhangHomodenZk::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_KZhangHomodenZk::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());

	sint size=0;
	// Количество выборок
	sint samples_number = (*SampVec).size();
	for (sint i=0; i < samples_number; i++)
	   size += (*SampVec)[i]->n();

	vector <sint> j_i(samples_number, 0);

	slongfloat  Fk, Fik;
	slongfloat Zk = 0.0,sum=0.0;

	vector<vector<sint> > R(samples_number);
	vector <elem> x(0);
	//вычислим ранги Rij

	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
			for (sint j = 0; j < ni; j++)
		{
			elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
		#pragma omp critical
		   {
			x.push_back(el);
		   }
		}
	}
   /*	x[0].y=2.0;
	x[1].y=5.0;
	x[2].y=6.0;
	x[3].y=7.0;
	x[4].y=1.0;
	x[5].y=3.0;
	x[6].y=4.0;
	x[7].y=8.0;    */

	sort( x.begin(), x.end(), ssort);
	// получаем ранг каждого наблюдения
	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
        vector < sfloat> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
			R[i][j]=buf[j] ;
	}

	for(sint k = 1; k <= size; k++)
	{
		sum=0.0;
		for(sint i = 0; i < samples_number; i++)
		{
			if(k==R[i][j_i[i]])
			{
				j_i[i]++;
				Fik = (j_i[i]-0.5)/sfloat((*SampVec)[i]->n());
			}
			else
				Fik = j_i[i]/sfloat((*SampVec)[i]->n());
			Fk = (k-0.5)/sfloat(size);
			if(  fabs(Fik/Fk)<1e-12)
				sum += ((*SampVec)[i]->n())*( (1-Fik)*log((1-Fik)/(1-Fk)));
			 else
				if(  fabs(1-Fik)<1e-12)
						sum += ((*SampVec)[i]->n())*(Fik * log(Fik/Fk)); 
				else
					sum += ((*SampVec)[i]->n())*(Fik * log(Fik/Fk) + (1-Fik)*log((1-Fik)/(1-Fk)));
		}

		Zk = dmax(Zk, sum);
	}
	   return Zk;
   }

   CriterionBase* Criterion_KZhangHomodenZk::Copy()
   {
	  return (new Criterion_KZhangHomodenZk());
   }

   Distribution * Criterion_KZhangHomodenZk::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
  }

 //----------------------------------------------------------------------------------------
 //---------------------------------------------------------------------------
 /// K-выборочный Критерий "Za"
  Criterion_KZhangHomodenZa::Criterion_KZhangHomodenZa()
   {
		_info.SetName("Критерий Жанга Za (К выборок)");
		_info.SetFileName("Критерий Жанга Za (К выборок)");
		_info.SetNote("Критерий Жанга Za (К выборок)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_KZhangHomodenZa::~Criterion_KZhangHomodenZa(){}

   void Criterion_KZhangHomodenZa::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_KZhangHomodenZa::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());

	sint size=0;
	// Количество выборок
	sint samples_number = (*SampVec).size();
	for (int i=0; i < samples_number; i++)
	   size += (*SampVec)[i]->n();

	vector <sint> j_i(samples_number, 0);
	slongfloat Fik;
	slongfloat Za = 0.0;

	vector<vector<sint> > R(samples_number);
	vector <elem> x(0);
	//вычислим ранги Rij

	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
			for (sint j = 0; j < ni; j++)
		{
			elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
		#pragma omp critical
		   {
			x.push_back(el);
		   }
		}
	}
	sort( x.begin(), x.end(), ssort);
// получаем ранг каждого наблюдения
	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
		vector < sfloat> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
					R[i][j]=buf[j] ;
	}

	for(sint k = 1; k <= size; k++)
		for(sint i = 0; i < samples_number; i++)
		{
			if(k==R[i][j_i[i]])
			{
				j_i[i]++;
				Fik = (j_i[i]-0.5)/sfloat((*SampVec)[i]->n());
			}
			else
				Fik = j_i[i]/sfloat((*SampVec)[i]->n());

			if(  fabs(Fik)<1e-12 || Fik==0)
				Za += ((*SampVec)[i]->n())*(1-Fik)*log(1-Fik)/((k-0.5)*(size-k+0.5));
			 else
				if(  fabs(1-Fik)<1e-12 || Fik==1 )
					Za += ((*SampVec)[i]->n()) * Fik * log(Fik)/((k-0.5)*(size-k+0.5));
				else
					Za += ((*SampVec)[i]->n()) *(Fik * log(Fik) + (1-Fik)*log(1-Fik) )/sfloat((k-0.5)*(size-k+0.5));
		}
	return -Za;
   }

   CriterionBase* Criterion_KZhangHomodenZa::Copy()
   {
	  return (new Criterion_KZhangHomodenZa());
   }

   Distribution * Criterion_KZhangHomodenZa::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
 }
   //----------------------------------------------------------------------------------------
//// Kвыборочный Критерий "Zc"
 Criterion_KZhangHomodenZc::Criterion_KZhangHomodenZc()
   {
		_info.SetName("Критерий Жанга Zc (К выборок)");
		_info.SetFileName("Критерий Жанга Zc (К выборок)");
		_info.SetNote("Критерий Жанга Zc (К выборок)");
		_info.SetGroup(CRITERIA_GROUP_HOMOGENEITY);
		_info.SetType(Nonparametric);
		FillRequiredParams();
   }

   Criterion_KZhangHomodenZc::~Criterion_KZhangHomodenZc(){}

   void Criterion_KZhangHomodenZc::FillRequiredParams()
   {
	  requiredParams.updateParameter(new ParameterDouble(ParameterName::SamplesList,2.0));
   }

   sfloat Criterion_KZhangHomodenZc::operator()(ParametersSet & params)
	{
	   std::vector<Sample*> * SampVec;
	   SampVec=(params.getParameter(ParameterName::SamplesList)->getVectorSample());

	sint size=0;
	// Количество выборок
	sint samples_number = (*SampVec).size();
	for (sint i=0; i < samples_number; i++)
	   size += (*SampVec)[i]->n();

	slongfloat Zc = 0.0;
	vector<vector<sint> > R(samples_number);
	vector <elem> x(0);

		//вычислим ранги Rij

	for(sint i = 0; i < samples_number; i++)
	{
		sint ni= (*SampVec)[i]->n();
		R[i].reserve(ni);
			for (sint j = 0; j < ni; j++)
		{
			elem el;
			el.ind=j;
			el.nSapmle=i;
			el.y = (*SampVec)[i]->xa(j);
			#pragma omp critical
		   {
			x.push_back(el);
		   }
		}
	}
	   /*	x[0].y=2.0;
	x[1].y=4.0;
	x[2].y=18.0;
	x[3].y=9.0;
	x[4].y=10.0;
	x[5].y=25.0;  */
	sort( x.begin(), x.end(), ssort);

   // получаем ранг каждого наблюдения
	for(sint i = 0; i < samples_number; i++)
	{
	sint ni= (*SampVec)[i]->n();
		vector < sfloat> buf(ni);
		for (sint j = 0; j < ni; j++)
			for (sint ii = 0; ii < size; ii++)
				 if (x[ii].nSapmle ==i && x[ii].ind==j)
					   buf[j]=ii+1;

		sort( buf.begin(), buf.end());
		R[i].reserve(ni);
		//R.push_back(buf);
		for (sint j = 0; j < ni; j++)
					R[i][j]=buf[j] ;
		buf.clear();
	}

	for(sint i = 0; i < samples_number; i++)
	{
		sfloat ni=(*SampVec)[i]->n();
		for(sint j = 0; j < ni; j++)
			Zc +=  log(ni/(j+1-0.5)-1.)*log(sfloat(size)/(R[i][j]-0.5)-1.);

	}
	return Zc/sfloat(size);
   }

   CriterionBase* Criterion_KZhangHomodenZc::Copy()
   {
	  return (new Criterion_KZhangHomodenZc());
   }

   Distribution * Criterion_KZhangHomodenZc::GetStatisticDistribution( ParametersSet * )
   {
		return NULL;
  }
   //----------------------------------------------------------------------------------------

   namespace{  	 RegTest r1(new Criterion_ZhangHomodenZk);
				 RegTest r2(new Criterion_ZhangHomodenZa);
			   	 RegTest r3(new Criterion_ZhangHomodenZc);
				 RegTest r4(new Criterion_KZhangHomodenZk);
				 RegTest r5(new Criterion_KZhangHomodenZa);
				 RegTest r6(new Criterion_KZhangHomodenZc);
			}



