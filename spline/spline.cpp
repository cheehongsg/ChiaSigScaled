#include "spline.h"
#include <iostream>
// ---------------------------------------------------------------------
// implementation part, which should be separated into a cpp file
// ---------------------------------------------------------------------


// band_matrix implementation
// -------------------------

band_matrix::band_matrix(long dim, long n_u, long n_l) {
	resize(dim, n_u, n_l);
}
void band_matrix::resize(long dim, long n_u, long n_l) {
	assert(dim>0);
	assert(n_u>=0);
	assert(n_l>=0);
	m_upper.resize(n_u+1);
	m_lower.resize(n_l+1);
	for(size_t i=0; i<m_upper.size(); i++) {
		m_upper[i].resize(dim);
	}
	for(size_t i=0; i<m_lower.size(); i++) {
		m_lower[i].resize(dim);
	}
}
long band_matrix::dim() const {
	if(m_upper.size()>0) {
		return m_upper[0].size();
	} else {
		return 0;
	}
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (long i, long j) {
	long k=j-i;       // what band is the entry
	assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
	assert( (-num_lower()<=k) && (k<=num_upper()) );
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if(k>=0)   return m_upper[k][i];
	else	    return m_lower[-k][i];
}
double band_matrix::operator () (long i, long j) const {
	long k=j-i;       // what band is the entry
	assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
	assert( (-num_lower()<=k) && (k<=num_upper()) );
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if(k>=0)   return m_upper[k][i];
	else	    return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double band_matrix::saved_diag(long i) const {
	assert( (i>=0) && (i<dim()) );
	return m_lower[0][i];
}
double & band_matrix::saved_diag(long i) {
	assert( (i>=0) && (i<dim()) );
	return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose() {
	long i_max,j_max;
	long j_min;
	double x;
	
	// preconditioning
	// normalize column i so that a_ii=1
	for(long i=0; i<this->dim(); i++) {
		assert(this->operator()(i,i)!=0.0);
		this->saved_diag(i)=1.0/this->operator()(i,i);
		j_min=std::max(0L,i-this->num_lower());
		j_max=std::min(this->dim()-1,i+this->num_upper());
		for(long j=j_min; j<=j_max; j++) {
			this->operator()(i,j) *= this->saved_diag(i);
		}
		this->operator()(i,i)=1.0;          // prevents rounding errors
	}
	
	// Gauss LR-Decomposition
	for(int k=0; k<this->dim(); k++) {
		i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
		for(int i=k+1; i<=i_max; i++) {
			assert(this->operator()(k,k)!=0.0);
			x=-this->operator()(i,k)/this->operator()(k,k);
			this->operator()(i,k)=-x;                         // assembly part of L
			j_max=std::min(this->dim()-1,k+this->num_upper());
			for(int j=k+1; j<=j_max; j++) {
				// assembly part of R
				this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
			}
		}
	}
}
// solves Ly=b
std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const {
	assert( this->dim()==(int)b.size() );
	std::vector<double> x(this->dim());
	long j_start;
	double sum;
	for(long i=0; i<this->dim(); i++) {
		sum=0;
		j_start=std::max(0L,i-this->num_lower());
		for(long j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
		x[i]=(b[i]*this->saved_diag(i)) - sum;
	}
	return x;
}
// solves Rx=y
std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const {
	assert( this->dim()==(int)b.size() );
	std::vector<double> x(this->dim());
	long j_stop;
	double sum;
	for(long i=this->dim()-1; i>=0; i--) {
		sum=0;
		j_stop=std::min(this->dim()-1,i+this->num_upper());
		for(long j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
		x[i]=( b[i] - sum ) / this->operator()(i,i);
	}
	return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
										  bool is_lu_decomposed) {
	assert( this->dim()==(int)b.size() );
	std::vector<double>  x,y;
	if(is_lu_decomposed==false) {
		this->lu_decompose();
	}
	y=this->l_solve(b);
	x=this->r_solve(y);
	return x;
}





// spline implementation
// -----------------------

spline::spline()
{
	
}

spline::~spline()
{
}

void spline::set_points(const std::vector<double>& x,
						const std::vector<double>& y, bool cubic_spline) {
	assert(x.size()==y.size());
	m_x=x;
	m_y=y;
	size_t n=x.size();
	// TODO sort x and y, rather than returning an error
	for(int i=0; i<n-1; i++) {
		assert(m_x[i]<m_x[i+1]);
	}
	
	if(cubic_spline==true) { // cubic spline interpolation
		// setting up the matrix and right hand side of the equation system
		// for the parameters b[]
		band_matrix A(n,1,1);
		std::vector<double>  rhs(n);
		for(int i=1; i<n-1; i++) {
			A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
			A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
			A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
			rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		}
		// boundary conditions, zero curvature b[0]=b[n-1]=0
		A(0,0)=2.0;
		A(0,1)=0.0;
		rhs[0]=0.0;
		A(n-1,n-1)=2.0;
		A(n-1,n-2)=0.0;
		rhs[n-1]=0.0;
		
		// solve the equation system to obtain the parameters b[]
		m_b=A.lu_solve(rhs);
		
		// calculate parameters a[] and c[] based on b[]
		m_a.resize(n);
		m_c.resize(n);
		for(int i=0; i<n-1; i++) {
			m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
			m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
			- 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
		}
	} else { // linear interpolation
		m_a.resize(n);
		m_b.resize(n);
		m_c.resize(n);
		for(int i=0; i<n-1; i++) {
			m_a[i]=0.0;
			m_b[i]=0.0;
			m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
		}
	}
	
	// for the right boundary we define
	// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
	double h=x[n-1]-x[n-2];
	// m_b[n-1] is determined by the boundary condition
	m_a[n-1]=0.0;
	m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
}

splineType spline::operator() (double x) const {
	size_t n=m_x.size();
	// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
	std::vector<double>::const_iterator it;
	it=std::lower_bound(m_x.begin(),m_x.end(),x);
	int idx=std::max( int(it-m_x.begin())-1, 0);
	
	double h=x-m_x[idx];
	double interpol;
	if(x<m_x[0]) {
		// extrapolation to the left
		interpol=((m_b[0])*h + m_c[0])*h + m_y[0];
	} else if(x>m_x[n-1]) {
		// extrapolation to the right
		interpol=((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];
	} else {
		// interpolation
		interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
	}
	return interpol;
}

splineType spline::cachedAt(int x) const {
	size_t n=m_x.size();
	// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
	std::vector<double>::const_iterator it;
	it=std::lower_bound(m_x.begin(),m_x.end(),(double)x);
	int idx=std::max( int(it-m_x.begin())-1, 0);
	
	double h=x-m_x[idx];
	double interpol;
	if(x<m_x[0]) {
		// extrapolation to the left
		interpol=((m_b[0])*h + m_c[0])*h + m_y[0];
	} else if(x>m_x[n-1]) {
		// extrapolation to the right
		interpol=((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];
	} else {
		// interpolation
		interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
	}
	return interpol;
}

