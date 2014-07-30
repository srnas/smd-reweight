/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 Giovanni Bussi

   smd-reweight is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <assert.h>
#include <cstdlib>

// Tool to allow 2D arrays
template<class T>
class Array2d{
public:
        Array2d(int row, int col):m_row(row),m_col(col),
        m_data((row!=0&&col!=0)?new T[row*col]:NULL){}
        Array2d(int row, int col,T val):m_row(row),m_col(col),
        m_data((row!=0&&col!=0)?new T[row*col]:NULL){
                for(int r=0;r<m_row;++r)for(int c=0;c<m_col;++c)(*this)[r][c]=val;
        }
        Array2d(const Array2d&src):m_row(src.m_row),m_col(src.m_col),
        m_data((src.m_row!=0&&src.m_col!=0)?new T[src.m_row*src.m_col]:NULL){
                for(int r=0;r<m_row;++r)for(int c=0;c<m_col;++c) (*this)[r][c] =
                        src[r][c];
        }
        ~Array2d(){if(m_data) delete []m_data;}
        inline T* operator[](int i) {return (m_data + (m_col*i));}
        inline T const* operator[](int i) const {return (m_data + (m_col*i));}
        inline int row()const{return m_row;}
        inline int col()const{return m_col;}
        inline void swap(Array2d&x){
          assert(row()==x.row() && col()==x.col());
          T*tmp; tmp=x.m_data; x.m_data=m_data; m_data=tmp;
        }
protected:
        Array2d& operator=(const Array2d&);
private:
        const int m_row;
        const int m_col;
        T* m_data;
};

using namespace std;

int main(int argc,char*argv[]){
// Default parameters
  double kT=2.49;
  int maxiter=10;
  int nbins=100;
  double tolerance=1e-4;

  double hmin,hmax;
  bool hminset=false;
  bool hmaxset=false;

  bool intbin=false;

  bool align_to_zero=false;

  string next="";
  for(int i=1;i<argc;i++){
    string a(argv[i]);
    istringstream is(a);
    if(next.length()>0){
      if     (next=="--nbins"){
        if(a!="int"){
          is>>nbins;
        } else intbin=true;
      }
      else if(next=="--maxiter")   is>>maxiter;
      else if(next=="--kt")        is>>kT;
      else if(next=="--tolerance") is>>tolerance;
      else if(next=="--hmin"){
        if(a!="auto"){
          is>>hmin;
          hminset=true;
        }
      }
      else if(next=="--hmax"){
        if(a!="auto"){
          is>>hmax;
          hmaxset=true;
        }
      }
      else assert(0);
      next="";
      continue;
    }
    if(a=="--nbins" || a=="--maxiter" || a=="--kt" || a=="--tolerance" || a=="--hmin" || a=="--hmax"){
      next=a;
    } else if(a=="--help" || a=="-h"){
      cout<<"\nUsage:\n"
          <<"  smd-reweight [flags] < file_with_concatenated_trajectories\n"
          <<"List of possible flags:\n"
          <<"--nbins     (default=100)    : number of bins in the analyzed CV ('int' means integer bins)\n"
          <<"--hmin      (default='auto') : minimum of the histogram ('auto' means min value)\n"
          <<"--hmax      (default='auto') : maximum of the histogram ('auto' means max value)\n"
          <<"--kt        (default=2.49)   : kt in energy units\n"
          <<"--maxiter   (default=10)     : maximum number of iterations in self-consistent cycle\n"
          <<"--tolerance (default=1e-4)   : tolerance in self-consistent cycle\n"
          <<"--align-to-zero              : align free-energy minumum\n\n"
          <<"Use '--nbins int' to have integer bins.\n"
          <<"The program expects from standard input a file with these columns:\n\n"
          <<"time pulled_cv position_of_restraint stiffness_of_restraint work analyzed_cv\n\n"
          <<"Lines beginning with '#' are ignored. First column is used to detect new trajectories,\n"
          <<"so that one can just concatenate them:\n\n"
          <<"cat COLVAR* | smd-reweight\n\n"
          <<"A suitable COLVAR file can be produced with this sample PLUMED 2.1 input:\n"
          <<"\n"
          <<"# pulled distance:\n"
          <<"d: DISTANCE ATOMS=1,2\n"
          <<"# moving restraint:\n"
          <<"r: MOVINGRESTRAINT AT0=0 AT1=1 STEP0=0 STEP1=100000 KAPPA0=20000.0\n"
          <<"# analyzed distance:\n"
          <<"a: DISTANCE ATOMS=3,4\n"
          <<"# print colvar file:\n"
          <<"PRINT ARG=d,r.d_cntr,r.d_kappa,r.d_work,a FILE=COLVAR\n"
          <<"\n"
          <<"All units should be consistent.\n\n";
      exit(0);
    } else if(a=="--align-to-zero"){
      align_to_zero=true;
    } else {
      cerr<<"ERROR: Unknown option "<<a<<"\n";
      exit(1);
    }
  }
  if(next.length()>0){
    cerr<<"ERROR: Needs argument for "<<next<<"\n";
    exit(1);
  }

  cout<<"# tolerance "<<tolerance<<endl;
  cout<<"# kt "<<kT<<endl;
  cout<<"# maximum number of iterations "<<maxiter<<endl;

  const double invkT=1.0/kT;
// number of frames
  int nframe=0;

// number of trajectories
  int ntraj=0;

// vector containing the entire input file
// stored line by line
  vector<string> file;
  {
    string line;
    double f,ff;
// read whole file:
    while(getline(cin,line)){
// skip comments and empty lines
      if(line.length()==0) continue;
      if(line.length()>0 && line[0]=='#') continue;
      istringstream is(line.c_str());
      is>>ff;
// first column is expected to be time and should decrease when a new trajectory is found
      if(nframe==0 && file.size()>0 && ff<f) nframe=file.size();
      f=ff;
      file.push_back(line);
    }
    ntraj=file.size()/nframe;
    cout<<"# number of frames "<<nframe<<endl;
    cout<<"# number of trajectories "<<ntraj<<endl;
// consistency check 
    if(file.size()%nframe!=0){
      cerr<<"ERROR: inconsistent trajectory files ("<<file.size()<<" lines for "<<nframe<<" frames\n";
      exit(1);
    }
  }
  
  Array2d<double> dist(nframe,ntraj);
  Array2d<double> w(nframe,ntraj);
  Array2d<double> dum(nframe,ntraj);
  Array2d<double> kappa(nframe,ntraj);
  Array2d<double> weights(nframe,ntraj);
  Array2d<double> previousWeights(nframe,ntraj);
  Array2d<double> ana(nframe,ntraj);

  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) 
  {
// read columns:
// time (ignored at this stage)
// position of restraint
// stiffness of restraint
// work
// variable to be analyzed
    istringstream is(file[itraj*nframe+iframe].c_str());
    double time;
    is>>time>>dist[iframe][itraj]>>dum[iframe][itraj]>>kappa[iframe][itraj]>>w[iframe][itraj]>>ana[iframe][itraj];
  }
  file.resize(0);


// find min and max for analyzed variable
  double ana_max=ana[0][0];
  double ana_min=ana[0][0];
  for(int iframe=0;iframe<nframe;iframe++) for(int itraj=0;itraj<ntraj;itraj++){
    if(ana[iframe][itraj]>ana_max){
      ana_max=ana[iframe][itraj];
    }
    if(ana[iframe][itraj]<ana_min){
      ana_min=ana[iframe][itraj];
    }
  }

  if(hminset) ana_min=hmin;
  if(hmaxset) ana_max=hmax;

  cout<<"# Histogram range "<<ana_min<<" "<<ana_max<<"\n";
  if(intbin){
    nbins=int(ana_max+0.5)-int(ana_min+0.5)+1;
  }
  cout<<"# Number of bins "<<nbins<<endl;
  vector<double> histo(nbins);

// Jarzynski calculation
  vector<double> jarz(nframe);
  for(int iframe=0;iframe<nframe;iframe++){
    double j=0.0;
    for(int itraj=0;itraj<ntraj;itraj++) j+=exp(-w[iframe][itraj]*invkT);
    jarz[iframe]=-kT*log(j);
  }

// nonequilibrium weights
  double noneq[nframe][ntraj];
  for(int iframe=0;iframe<nframe;iframe++) for(int itraj=0;itraj<ntraj;itraj++)
    noneq[iframe][itraj]=exp(-w[iframe][itraj]*invkT);


// non-self-consistent calculation (just as a check)
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
      double e=dist[iframe][itraj]-dum[iframe][itraj];
      weights[iframe][itraj]=exp(-(w[iframe][itraj]-0.5*kappa[iframe][itraj]*e*e)*invkT);
  }

// normalization of weights
  double norm=0.0;
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) norm+=weights[iframe][itraj];
  double invnorm=1.0/norm;
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) weights[iframe][itraj]*=invnorm;

// plain jarzinsky used for first guess of F
  vector<double> F=jarz;
  {
    ofstream of("F.jarz");
    for(int iframe=0;iframe<nframe;iframe++) of<<dum[iframe][0]<<" "<<F[iframe]<<endl;
  }

// self-consistent iterations
for(int iter=0;iter<maxiter;iter++){

// store past weights
  previousWeights.swap(weights);

// compute weights
  double norm=0.0;
#pragma omp parallel for
  for(int iframe=0;iframe<nframe;iframe++){
    vector<double> tmp(nframe);
    for(int itraj=0;itraj<ntraj;itraj++){
      double j=0;
      double kk=-0.5*kappa[iframe][itraj];
      double dd=dist[iframe][itraj];
      for(int jframe=0;jframe<nframe;jframe++){
        double e=dd-dum[jframe][itraj];
        tmp[jframe]=kk*e*e*invkT+F[jframe]*invkT;
      }
      for(int jframe=0;jframe<nframe;jframe++) j+=exp(tmp[jframe]);
#pragma omp critical
      {
        weights[iframe][itraj]=noneq[iframe][itraj]*exp(F[iframe]*invkT)/j;
        norm+=weights[iframe][itraj];
      }
    }
  }
  double invnorm=1.0/norm;
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) weights[iframe][itraj]*=invnorm;

// recompute weighted histogram
  histo.assign(histo.size(),0.0);
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
    int i=(int) (0.5+(ana[iframe][itraj]-ana_min)/(ana_max-ana_min)*(histo.size()-1));
    if(i==histo.size())i=histo.size()-1;
    assert(i>=0 && i<histo.size());
    histo[i]+=weights[iframe][itraj];
  }
  if(align_to_zero){
    double zero=histo[0];
    for(int i=0;i<histo.size();i++) if(histo[i]>zero)zero=histo[i];
    for(int i=0;i<histo.size();i++) histo[i]/=zero;
  }
  {
// this is a small file, so we rewrite it at every iteration
    ofstream of("histo");
    for(int i=0;i<histo.size();i++) of<<ana_min+i*(ana_max-ana_min)/(histo.size()-1)<<" "
                                      <<-kT*log(histo[i])<<endl;
  }

// update estimated F
#pragma omp parallel for
  for(int jframe=0;jframe<nframe;jframe++){
    vector<double> tmp(nframe);
    double j=0;
    for(int itraj=0;itraj<ntraj;itraj++){
      for(int iframe=0;iframe<nframe;iframe++) {
        double e=dist[iframe][itraj]-dum[jframe][itraj];
        tmp[iframe]=-0.5*kappa[iframe][itraj]*e*e*invkT;
      }
      for(int iframe=0;iframe<nframe;iframe++) j+=weights[iframe][itraj]*exp(tmp[iframe]);
    }
#pragma omp critical
    F[jframe]=-kT*log(j);
  }

// compute error
  double eps=0.0;
#pragma omp parallel for
  for(int iframe=0;iframe<nframe;iframe++){
    double tmp=0;
    for(int itraj=0;itraj<ntraj;itraj++)
    {
      double e=kT*log(weights[iframe][itraj]/previousWeights[iframe][itraj]);
      tmp+=e*e;
    }
#pragma omp critical
    eps+=tmp;
  }
  eps/=(nframe*ntraj);

  cout<<"# iteration "<<iter<<" error "<<sqrt(eps)<<endl;

  if(sqrt(eps)<tolerance) break;

// final F
  {
    ofstream of("F.final");
    for(int iframe=0;iframe<nframe;iframe++) of<<dum[iframe][0]<<" "<<F[iframe]<<endl;
  }
// final -kTlogweight
  {
    ofstream of("weight.final");
    for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) of<<-kT*log(weights[iframe][itraj])<<endl;
  }


}
  
  return 0;
}
