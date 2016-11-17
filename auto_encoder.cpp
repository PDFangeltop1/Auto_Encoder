#include <iostream>
#include <ctime>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#define PI 3.141592653589

using namespace std;
double rand(double min, double max){
  return min+(max-min)*rand()/(RAND_MAX+1.0);
}
double normal(double x, double miu, double sigma){
  return exp(-1.0*(x-miu)*(x-miu)/(2*sigma*sigma))/(sqrt(2*PI)*sigma);
}
double randn(double miu, double sigma, double min, double max){
  double x,y, dScope;
  do{
    x=rand(min,max);
    y=normal(x,miu,sigma);
    dScope=rand(0.0,normal(miu,miu,sigma));
  }while(dScope>y);
  return x;
}

void initMat(vector<vector<double> >& mat,int x,int y,string kind){
  mat.clear();
  for(int i=0; i<x; i++){
    vector<double> tmp;
    for(int j=0; j<y; j++){
      if(kind == "normal"){
	tmp.push_back(randn(0,0.5,-0.2,0.2));
      }else if(kind == "zero"){
	tmp.push_back(0.0);
      }
    }
    mat.push_back(tmp);
  }
}

class Matrix{
public:
  vector<vector<double> >mat;

  //double a[][];
  int h;
  int w;
  Matrix(){}
  Matrix(int h, int w,string kind){
    this->h = h;
    this->w = w;
    initMat(this->mat,h,w,kind);
  }
  Matrix(const Matrix& x){
    this->h = x.h;
    this->w = x.w;
    initMat(this->mat,x.h,x.w,"zero");
    for(int i=0; i<x.h; i++){
      for(int j=0; j<x.w; j++){
	mat[i][j] = x.mat[i][j];
      }
    }
  }
  Matrix(const vector<double>& x){
    this->h = 1;
    this->w = x.size();
    initMat(this->mat,1,x.size(),"zero");
    for(int i=0; i<x.size(); i++){
      mat[0][i] = x[i];
    }
  }
  double norm(){
    double ans = 0;
    for(int i=0; i<mat.size(); i++){
      for(int j=0; j<mat[0].size(); j++){
	ans += mat[i][j]*mat[i][j];
      }
    }
    return sqrt(ans);
  }
  Matrix& operator* (const Matrix& y){
    assert(h == y.h);
    assert(w == y.w);
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] *= y.mat[i][j];
      }
    }
    return *this;
  }
  Matrix& operator+ (const Matrix& y){
    assert(h == y.h);
    assert(w == y.w);
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] += y.mat[i][j];
      }
    }
    return *this;
  }
  Matrix& operator- (const Matrix& y){
    assert(h == y.h);
    assert(w == y.w);
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] -= y.mat[i][j];
      }
    }
    return *this;
  }
  Matrix& operator/ (const Matrix& y){
    assert(h == y.h);
    assert(w == y.w);
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] /= y.mat[i][j];
      }
    }
    return *this;
  }
  Matrix& mat_sqrt(){
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] = sqrt(this->mat[i][j]);
      }
    }
    return *this;
  }
  Matrix& scalar_add(double a){
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] += a;
      }
    }
    return *this;    
  }
  Matrix& scalar_mut(double a){
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	this->mat[i][j] *= a;
      }
    }
    return *this;    
  }
  double inter(const Matrix& y){
    assert(h == y.h);
    assert(w == y.w);
    double ans = 0;
    for(int i=0; i<h; i++){
      for(int j=0; j<w ;j++){
	ans += this->mat[i][j]*y.mat[i][j];
      }
    }
    return ans;
  }
  Matrix* transpose(){
    Matrix* ans = new Matrix(w,h,"zero");
    for(int i=0; i<w; i++){
      for(int j=0; j<h; j++){
	ans->mat[i][j] = this->mat[j][i];
      }
    }
    return ans;
  }
  Matrix outer(Matrix& y){
    assert(h == 1);
    assert(y.h == 1);
    Matrix ans = Matrix(w,y.w, "zero");
    for(int i=0 ;i<w; i++){
      for(int j=0; j<y.w; j++){
	ans.mat[i][j] = mat[0][i]*y.mat[0][j];
      }
    }
    return ans;
  }
  Matrix dot(Matrix& W){
    //W: D*M
    //x: N*D
    assert(w == W.h);
    Matrix ans = Matrix(h,W.w,"zero");
    for(int i=0; i<h; i++){
      for(int j=0; j<W.w; j++){
	for(int k=0; k<w; k++){
	  ans.mat[i][j] += mat[i][k]*W.mat[k][j];
	}
      }
    }
    return ans;
  }
  
};

void parameter_to_file(ofstream& myfile, const Matrix& mat){
  for(int i=0; i<mat.h; i++){
    for(int j=0; j<mat.w; j++){
      myfile<<mat.mat[i][j]<<" ";
    }
    myfile<<"\n";
  }
}

double min_square_loss(const Matrix& x, const Matrix& y){
  assert(x.h == 1);
  assert(y.h == 1);
  assert(x.w == y.w);
  double ans = 0;
  for(int i=0; i<x.w; i++){
    ans += (x.mat[0][i]-y.mat[0][i])*(x.mat[0][i]-y.mat[0][i]);
  }
  return ans*0.5;
}

double abss(double x){
  if(x< 0) x = -1.0*x;
  return x;
}

double maxx(double x, double y){
  if(x > y){
    return x;
  }else{
    return y;
  }
}

double rel_error_mat(Matrix& x, Matrix& y){
  assert(x.h == y.h);
  assert(x.w == y.w);
  double ans = 0;
  for(int i=0; i<x.h; i++){
    for(int j=0; j<x.w; j++){
      double tmp = abss(x.mat[i][j]-y.mat[i][j])/maxx(1e-8,abss(x.mat[i][j])+abss(y.mat[i][j]));
      if(tmp > ans){
	ans = tmp;
      }
    }
  }
  return ans;
}

Matrix affine_forward(Matrix& x, Matrix& W, Matrix& b){
  //x: N*D
  //W: D*M
  //b: M
  assert(x.w == W.h);
  assert(W.w == b.w);
  Matrix ans = x.dot(W)+b;
  return ans;
}

void affine_backward(Matrix& dout, 
		     Matrix& x, Matrix& W, Matrix& b,
		     Matrix& dx, Matrix& dW, Matrix& db){
  
  Matrix* tran = W.transpose();
  dx = dout.dot(*tran);
  dW = dW + x.outer(dout);
  db = db + dout;
  delete tran;
}

void sgd_update_parameter(Matrix& W, Matrix& dW, 
			     double alpha){
  assert(W.h == dW.h);
  assert(W.w == dW.w);
  W = W - dW.scalar_mut(alpha);
}

class Adam{
public:
  Matrix* m;
  Matrix* v;
  double learning_rate;
  double beta1;
  double beta2;
  double epsilon;
  double t;
  Adam(double lr,double be1,double be2, double e, Matrix& x){
    this->learning_rate = lr;
    this->beta1 = be1;
    this->beta2 = be2;
    this->epsilon = e;
    this->t = 0;
    m = new Matrix(x.h, x.w, "zero");
    v = new Matrix(x.h, x.w, "zero");
  }
  ~Adam(){
    delete m;
    delete v;
  }
};
void Adam_update_parameter(Matrix& W,Matrix& dW, Adam& config){
  config.t += 1;
  (*config.m) = config.m->scalar_mut(config.beta1) + dW.scalar_mut(1-config.beta1);
  Matrix dW_2 = dW*dW;
  (*config.v) = config.v->scalar_mut(config.beta2) + dW_2.scalar_mut(1-config.beta2);
  Matrix x1 = config.m->scalar_mut(config.learning_rate);
  Matrix sqrt_tmp = config.v->mat_sqrt();
  Matrix x2 = sqrt_tmp.scalar_add(config.epsilon);
  W = W - x1/x2;
}

class MomentumSGD{
public:
  double learning_rate;
  double momentum;
  Matrix* velocity;
  MomentumSGD(double lr, double mu, Matrix& x){
    this->learning_rate = lr;
    this->momentum = mu;  
    velocity = new Matrix(x.h, x.w,"zero");
  }
  ~MomentumSGD(){
    delete velocity;
  }
};

void  MomentumSGD_update_parameter(Matrix& W, Matrix& dW,MomentumSGD& config){
  *config.velocity = config.velocity->scalar_mut(config.momentum) - dW.scalar_mut(config.learning_rate);
  W = W + *config.velocity;
}

class AutoEncoder{
public:
  int d,h;
  Matrix *W1;  Matrix *W2;
  Matrix *b1;  Matrix *b2;
  Matrix *dW1;  Matrix *dW2;
  Matrix *db1;  Matrix *db2;
    
  AutoEncoder(){}
  ~AutoEncoder(){
    //cout<<" xi gou han shu"<<endl;
    delete W1;    delete W2;    delete b1;
    delete b2;    delete dW1;    delete dW2;
    delete db1;    delete db2;
  }
  AutoEncoder(int d,int h);
  AutoEncoder& operator= (const AutoEncoder& x);
  double get_loss(Matrix& x);
  void update_parameter_momentum(double alpha,int iter, 				   
				 MomentumSGD& configW1, MomentumSGD& configW2,
				 MomentumSGD& configb1, MomentumSGD& configb2);
  void update_parameter_sgd(double alpha, int iter);
  void update_parameter_Adam(double alpha, Adam& configW1, Adam& configW2,
			     Adam& configb1, Adam& configb2);
  
  void grad_zeros();
  void gradient_check(Matrix& x);
  void eval_numerical_gradient_one_layer(Matrix& x,Matrix& W, Matrix&b,
					 Matrix& dmat_num, double delta, 
					 Matrix& dout, Matrix& change);
  void eval_numerical_gradient_two_layers(Matrix& x,Matrix& W, Matrix&b,
					  Matrix&W2, Matrix& b2,
					  Matrix& dmat_num, double delta, 
					  Matrix& dout, Matrix& change);

  void load_parameter(string file_name);
  void save_parameter(double loss,string x);
};

AutoEncoder& AutoEncoder::operator=(const AutoEncoder& x){
  d = x.d;
  h = x.h;
  *W1 = *x.W1;
  *W2 = *x.W2;
  *b1 = *x.b1;
  *b2 = *x.b2;
  return *this;
}

AutoEncoder::AutoEncoder(int d, int h){
  this->d = d;
  this->h = h;
  this->b1 = new Matrix(1,h,"zero");
  this->b2 = new Matrix(1,d,"zero");
  this->W1 = new Matrix(d,h,"normal");
  this->W2 = new Matrix(h,d,"normal");
  this->db1 = new Matrix(1,h,"zero");
  this->db2 = new Matrix(1,d,"zero");
  this->dW1 = new Matrix(d,h,"zero");
  this->dW2 = new Matrix(h,d,"zero");
}

void AutoEncoder::grad_zeros(){
  initMat(this->db1->mat,1,h,"zero");
  initMat(this->db2->mat,1,d,"zero");
  initMat(this->dW1->mat,d,h,"zero");
  initMat(this->dW2->mat,h,d,"zero");  
}

void AutoEncoder::load_parameter(string file_name){
  string line;
  const char *cfile_name = file_name.c_str();
  ifstream myfile;
  myfile.open(cfile_name);
  if(myfile.is_open()){
    int count = 0;
    while(getline(myfile,line)){
      stringstream ss(line);
      if(count < 5){
	for(int i=0; i<10; i++){
	  double x;
	  ss>>x;
	  this->W1->mat[i][count] = x;
	}
      }	else if(count>=5 && count < 15){
	for(int i=0; i<5; i++){
	  double x;
	  ss>>x;
	  this->W2->mat[i][count-5] = x;
	}
      }	else if(count == 15){
	for(int i=0; i<5; i++){
	  double x;
	  ss>>x;
	  this->b1->mat[0][i] = x;
	}	
      }	else if(count == 16){
	for(int i=0; i<10; i++){
	  double x;
	  ss>>x;
	  this->b2->mat[0][i] = x;
	}	
      }      
      count += 1;
    }
    myfile.close();
    cout<<"load parameter successfully"<<endl;
  }else{
    cout<<" unable to open file"<<endl;
  }  
}
void AutoEncoder::save_parameter(double loss,string dir_name){
  stringstream ss;
  ss<<loss;
  string file_name = "parameters_" + dir_name+ "/" + ss.str() + ".para";
  const char *cfile_name = file_name.c_str();
  ofstream myfile;
  myfile.open(cfile_name);
  if(myfile.is_open()){
    parameter_to_file(myfile,*this->W1);
    parameter_to_file(myfile,*this->W2);
    parameter_to_file(myfile,*this->b1);
    parameter_to_file(myfile,*this->b2);
  }else{
    cout<<"unable to open file"<<endl;
  }
}

void AutoEncoder::eval_numerical_gradient_one_layer(Matrix& x,Matrix& W, Matrix&b,
						   Matrix& dmat_num, double delta, Matrix& dout, Matrix& change){
  for(int i=0; i<dmat_num.mat.size(); i++){
    for(int j=0; j<dmat_num.mat[i].size(); j++){
      double old_value = change.mat[i][j];
      change.mat[i][j] = old_value + delta;
      Matrix fxp = affine_forward(x,W,b);
      change.mat[i][j] = old_value - delta;
      Matrix fxm = affine_forward(x,W,b);
      change.mat[i][j] = old_value;
      dmat_num.mat[i][j] = (fxp-fxm).inter(dout)/(2*delta);
    }
  }
}
void AutoEncoder::eval_numerical_gradient_two_layers(Matrix& x,Matrix& W, Matrix&b,
						     Matrix&W2, Matrix& b2,
						     Matrix& dmat_num, double delta, 
						     Matrix& dout, Matrix& change){
  for(int i=0; i<dmat_num.mat.size(); i++){
    for(int j=0; j<dmat_num.mat[i].size(); j++){
      double old_value = change.mat[i][j];
      change.mat[i][j] = old_value + delta;
      Matrix fxp = affine_forward(x,W,b);
      fxp = affine_forward(fxp,W2,b2);
      change.mat[i][j] = old_value - delta;
      Matrix fxm = affine_forward(x,W,b);
      fxm = affine_forward(fxm,W2,b2);
      change.mat[i][j] = old_value;
      dmat_num.mat[i][j] = (fxp-fxm).inter(dout)/(2*delta);
    }
  }
}

void AutoEncoder::gradient_check(Matrix& x){
  Matrix W = Matrix(10,5,"normal");
  Matrix b = Matrix(1,5,"zero");
  Matrix W2 = Matrix(5,10, "normal");
  Matrix b2 = Matrix(1,10,"zero");
  Matrix dW = Matrix(10,5,"zero");
  Matrix db = Matrix(1,5,"zero");
  Matrix dx = Matrix(1,10,"zero");
  Matrix dW2 = Matrix(5,10,"zero");
  Matrix db2 = Matrix(1,10,"zero");

  
  Matrix h = affine_forward(x,W,b);
  Matrix y = affine_forward(h,W2,b2);
  Matrix dy =y-x;
  Matrix dh = Matrix(1,5,"zero");
  affine_backward(dy,h,W2,b2,dh,dW2,db2);
  affine_backward(dh,x,W,b,dx,dW,db);
  
  Matrix dW_num = Matrix(10,5,"zero");
  Matrix db_num = Matrix(1,5,"zero");
  Matrix dx_num = Matrix(1,10,"zero");
  Matrix dW2_num = Matrix(5,10,"zero");
  Matrix db2_num = Matrix(1,10,"zero");
  Matrix dh_num = Matrix(1,5,"zero");;
  double delta=0.00001;
  eval_numerical_gradient_two_layers(x,W,b,W2,b2,db_num,delta,dy,b);
  eval_numerical_gradient_two_layers(x,W,b,W2,b2,dW_num,delta,dy,W);
  eval_numerical_gradient_two_layers(x,W,b,W2,b2,db2_num,delta,dy,b2);
  eval_numerical_gradient_two_layers(x,W,b,W2,b2,dW2_num,delta,dy,W2);

  //cout<<" error of x "<< rel_error_vec(dx_num, dx)<<endl;
  cout<<" error of b1 "<< rel_error_mat(db_num, db)<<endl;
  cout<<" error of W1 "<< rel_error_mat(dW_num, dW)<<endl;
  cout<<" error of b2 "<< rel_error_mat(db2_num, db2)<<endl;
  cout<<" error of W2 "<< rel_error_mat(dW2_num, dW2)<<endl;
}

double AutoEncoder::get_loss(Matrix& x){
  // Propagate
  Matrix h = affine_forward(x,*this->W1,*this->b1);
  Matrix y = affine_forward(h,*this->W2,*this->b2);
  double loss = min_square_loss(x,y);
  
  // Back propagate 
  y = y-x;
  Matrix grad_y = Matrix(y);
  Matrix grad_h;
  affine_backward(grad_y,h,*this->W2,*this->b2,grad_h,*this->dW2,*this->db2);
  Matrix grad_x;
  affine_backward(grad_h,x,*this->W1,*this->b1,grad_x,*this->dW1,*this->db1);
  return loss;
}

void AutoEncoder::update_parameter_sgd(double alpha, int iter){
  double dW2_norm = this->dW2->norm();
  double W2_norm = this->W2->norm();
  double tmp_dW2_norm = W2_norm;
  cout<< " W2 norm is: "<<W2_norm <<" dW2 norm is :"<< dW2_norm*alpha <<"  ";  
  sgd_update_parameter(*this->W2,*this->dW2,alpha);
  sgd_update_parameter(*this->W1,*this->dW1,alpha);
  sgd_update_parameter(*this->b2,*this->db2,alpha);
  sgd_update_parameter(*this->b1,*this->db1,alpha);
}

void AutoEncoder::update_parameter_Adam(double alpha,
					Adam& configW1, Adam& configW2,
					Adam& configb1, Adam& configb2){
  double dW2_norm = this->dW2->norm();
  double W2_norm = this->W2->norm();
  double tmp_dW2_norm = W2_norm;
  cout<< " W2 norm is: "<<W2_norm <<" dW2 norm is :"<< dW2_norm*alpha <<"  ";
  Adam_update_parameter(*this->W2,*this->dW2,configW2);
  Adam_update_parameter(*this->W1,*this->dW1,configW1);
  Adam_update_parameter(*this->b2,*this->db2,configb2);
  Adam_update_parameter(*this->b1,*this->db1,configb1);
}
void AutoEncoder::update_parameter_momentum(double alpha,int iter, 
					    MomentumSGD& configW1, MomentumSGD& configW2,
					    MomentumSGD& configb1, MomentumSGD& configb2){
  double dW2_norm = this->dW2->norm();
  double W2_norm = this->W2->norm();
  double tmp_dW2_norm = W2_norm;
  cout<< " W2 norm is: "<<W2_norm <<" dW2 norm is :"<< dW2_norm*alpha <<"  ";
  MomentumSGD_update_parameter(*this->W2,*this->dW2,configW2);
  MomentumSGD_update_parameter(*this->W1,*this->dW1,configW1);
  MomentumSGD_update_parameter(*this->b2,*this->db2,configb2);
  MomentumSGD_update_parameter(*this->b1,*this->db1,configb1);
}

class Solver{
public:
  double eta;
  int epoch;

  //parameter for momentum
  double mu;
  MomentumSGD *b1,*b2, *W1,*W2;
  //parameter for Adam
  double beta1, beta2, epsilon;
  Adam *aW1,*aW2, *ab1,*ab2;
  
  AutoEncoder* model;
  vector<vector<double> > training_data;
  Solver(AutoEncoder* encoder, double mu,double eta, double beta1,double beta2, double epsilon,int epoch);
  ~Solver();
  void load_training_data(string file_name);
  void train(int optimization_rule);
  void change_momentum(double new_mu);
  double get_training_data_loss();
};

Solver::Solver(AutoEncoder* encoder,double mu,double eta,double beta1,double beta2, double epsilon,int epoch){
  this->model = encoder;
  this->eta = eta;
  this->epoch = epoch;
  this->mu = mu;
  this->b1 = new MomentumSGD(eta,mu,*encoder->b1);   
  this->b2 = new MomentumSGD(eta,mu,*encoder->b2);   
  this->W1 = new MomentumSGD(eta,mu,*encoder->W1);   
  this->W2 = new MomentumSGD(eta,mu,*encoder->W2);   
  this->ab1 = new Adam(eta,beta1,beta2,epsilon,*encoder->b1);   
  this->ab2 = new Adam(eta,beta1,beta2,epsilon,*encoder->b2);   
  this->aW1 = new Adam(eta,beta1,beta2,epsilon,*encoder->W1);   
  this->aW2 = new Adam(eta,beta1,beta2,epsilon,*encoder->W2);   
}
Solver::~Solver(){
  delete this->model;
  delete this->b1;
  delete this->b2;
  delete this->W1;
  delete this->W2;
  delete this->ab1;
  delete this->ab2;
  delete this->aW1;
  delete this->aW2;
  
}

void Solver::change_momentum(double new_mu){
  this->b1->momentum = new_mu;
  this->b2->momentum = new_mu;
  this->W1->momentum = new_mu;
  this->W2->momentum = new_mu;
}
double Solver::get_training_data_loss(){
  double loss = 0;
  for(int j=0; j<this->training_data.size(); j++){
    Matrix tmp = Matrix(this->training_data[j]);
    loss += this->model->get_loss(tmp);
  }
  loss /= this->training_data.size();
  return loss;
}
void Solver::load_training_data(string file_name){
  string line;
  ifstream myfile("dataset.dat");
  if(myfile.is_open()){
    int count = 0;
    while(getline(myfile,line)){
      stringstream ss(line);
      if(count == 0){
	count += 1;
      }else{
	vector<double> tmp;
	for(int i=0; i<10; i++){
	  double x;
	  ss>>x;
	  tmp.push_back(x);
	}
	this->training_data.push_back(tmp);
      }
    }
    myfile.close();
    cout<<"load data successfully"<<endl;
  }else{
    cout<<" unable to open file"<<endl;
  }
}

void Solver::train(int optimization_rule){
  Matrix tmp = Matrix(this->training_data[0]);
  this->model->gradient_check(tmp);
  int flag = 1;
  double loss = 0;
  double prev_loss=0;
  for(int i=0; i<this->epoch; i++){
    this->model->grad_zeros();
    for(int j=0; j<this->training_data.size(); j++){
      Matrix tmp = Matrix(this->training_data[j]);
      loss += this->model->get_loss(tmp);
    }
    loss /= this->training_data.size();    
    if(optimization_rule == 1){
      //Adam  
      this->model->update_parameter_Adam(this->eta,*this->aW1,*this->aW2,*this->ab1,*this->ab2);
      cout<<" adam ";
    }else if(optimization_rule == 2){
      // Momentum_sgd 
      this->model->update_parameter_momentum(this->eta, i, *this->W1,*this->W2, *this->b1,*this->b2);
      if(loss < 0.12){
        this->change_momentum(0.9);
      }
      cout<<" momentum "<<endl;
    }else if(optimization_rule == 3){
      //sgd
      this->model->update_parameter_sgd(this->eta,i);
      cout<< "sgd"<<endl;
    }
    if(i >= 1 && i%1000 == 0 && abss(loss-prev_loss) > 1e-6){
      this->eta *= 0.95;
    }
    if(loss < prev_loss && loss < 0.111){
      this->model->save_parameter(loss,"charanda03");
    }
    cout<<" epoch : "<<i<<" ";
    cout<<" loss : "<<loss<<endl;
    prev_loss = loss;
    loss = 0;
  }
}

int main(){  
  srand((unsigned)time(NULL));
  AutoEncoder* t = new AutoEncoder(10,5);
  //t->load_parameter("chainer_param_today");
  Solver s = Solver(t,0.5,1e-3,0.9,0.999,1e-8,200);   //encoder, mu , eta, beta1,beta2,epsilon,epoch
  s.load_training_data("dataset.dat");
  s.train(1);
  //double loss = s.get_training_data_loss();
  //cout<<" loss is "<<loss<<endl;
  return 0;
}

