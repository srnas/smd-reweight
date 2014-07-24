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
  const double kT=2.49;
  const double invkT=1.0/kT;
  const int maxiter=10;
  const int nbins=100;
  const double tolerance=1e-4;

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
    assert(file.size()%nframe==0);
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

  vector<double> histo(nbins);

// find min and max for analyzed variable
  double ana_max=ana[0][0];
  double ana_min=ana[0][0];
  for(int iframe=0;iframe<nframe;iframe++) for(int itraj=0;itraj<ntraj;itraj++){
    if(ana[iframe][itraj]>ana_max){
      ana_max=ana[iframe][itraj]+1e-5;
    }
    if(ana[iframe][itraj]<ana_min){
      ana_min=ana[iframe][itraj]-1e-5;
    }
  }

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
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) weights[iframe][itraj]/=norm;

// change this if you want the non-consistent weighted histogram
  if(false){
    histo.assign(histo.size(),0.0);
    for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
      int i=(int) ((ana[iframe][itraj]-ana_min)/(ana_max-ana_min)*histo.size());
      if(i==histo.size())i=histo.size()-1;
      assert(i>=0 && i<histo.size());
      histo[i]+=weights[iframe][itraj];
    }

    ofstream of("histo-non-consistent");
    for(int i=0;i<histo.size();i++) of<<ana_min+(ana_max-ana_min)*(i+0.5)/histo.size()<<" "
                                      <<-kT*log(histo[i])<<endl;
    ofstream ofw("weights-non-consistent");
    for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
      ofw<<-kT*log(weights[iframe][itraj])<<endl;
    }
  }

// plain jarzinsky used for first guess of F
  vector<double> F=jarz;
  {
    ofstream of("F.jarz");
    for(int iframe=0;iframe<nframe;iframe++) of<<dum[iframe][0]<<" "<<F[iframe]<<endl;
  }

// self-consistent iterations
for(int iter=0;iter<maxiter;iter++){

// store past weights
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) previousWeights[iframe][itraj]=weights[iframe][itraj];

// compute weights
  double norm=0.0;
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) {
    double j=0;
    for(int jframe=0;jframe<nframe;jframe++){
      double e=dist[iframe][itraj]-dum[jframe][itraj];
      j+=exp(-0.5*kappa[iframe][itraj]*e*e*invkT+F[jframe]*invkT);
    }
    weights[iframe][itraj]=noneq[iframe][itraj]*exp(F[iframe]*invkT)/j;
    norm+=weights[iframe][itraj];
  }
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) weights[iframe][itraj]/=norm;

// recompute weighted histogram
  histo.assign(histo.size(),0.0);
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
    int i=(int) ((ana[iframe][itraj]-ana_min)/(ana_max-ana_min)*histo.size());
    if(i==histo.size())i=histo.size()-1;
    assert(i>=0 && i<histo.size());
    histo[i]+=weights[iframe][itraj];
  }
  {
// this is a small file, so we rewrite it at every iteration
    ofstream of("histo");
    for(int i=0;i<histo.size();i++) of<<ana_min+(ana_max-ana_min)*(i+0.5)/histo.size()<<" "
                                      <<-kT*log(histo[i])<<endl;
  }

// update estimated F
  for(int jframe=0;jframe<nframe;jframe++){
    double j=0;
    for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++) {
      double e=dist[iframe][itraj]-dum[jframe][itraj];
      j+=weights[iframe][itraj]*exp(-0.5*kappa[iframe][itraj]*e*e*invkT);
    }
    F[jframe]=-kT*log(j);
  }

// compute error
  double eps=0.0;
  for(int itraj=0;itraj<ntraj;itraj++) for(int iframe=0;iframe<nframe;iframe++){
    double e=kT*log(weights[iframe][itraj]/previousWeights[iframe][itraj]);
    eps+=e*e;
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
