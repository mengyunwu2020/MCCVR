#include <RcppArmadillo.h> 

using namespace Rcpp;
using namespace std;
using namespace arma;


//// [[Rcpp::export]]

arma::mat Newtonmultinomial(arma::mat GZ, arma::mat B, arma::mat alpha, arma::mat beta, double mu, arma::mat YL, double h, arma::mat Zini, arma::mat Zini1, arma::mat Zini2, arma::mat Zini3)
{
	int maxiter = 20;
	int n = GZ.n_rows;
	arma::mat Z = Zini;
	arma::mat Z1 = Zini1;
	arma::mat Z2 = Zini2;
	arma::mat Z3 = Zini3;
	arma::mat MM = -GZ / h + B*beta + ones(n, 1)*alpha;
	arma::mat expZ, dz, Hz;
	
	int iter = 0;
	for (iter = 0; iter < maxiter; iter++)
	{
		expZ = exp(Z);
		arma::mat sumexpZ = exp(Z)+ exp(Z1) + exp(Z2) + exp(Z3);
		dz = -mu*(YL - expZ / sumexpZ) / n + h*(Z - MM);
		Hz = mu*(expZ % (sumexpZ - expZ)) / square(sumexpZ) / n + h;
		Z = Z - dz / Hz;
	}
	return(Z);
}


//// [[Rcpp::export]]
// min_B {|Y-XB|^2 + lam*|B|}
arma::mat MGlasso_C(arma::mat Y, arma::mat X, arma::vec lam, arma::mat B0, double conv, int maxiter) {
  int  p=X.n_cols, iter=0, j;  // n=Y.n_rows, q=Y.n_cols,
  double  diff=10*conv, l2B1;
  rowvec sh; 
  arma::mat mat1=eye(p,p), B1, res1, res1j, XRj; 
  
  if (lam.size() == 1) {lam = as_scalar(lam)*ones(p);}
  sh = sum(square(X), 0);
  if (B0.is_finite()) {
    B1 = B0;
  } else {
    Rcpp::Rcout <<"B need to be initialized"<< std::endl;
    B1 = solve(X.t()*X, X.t()*Y, solve_opts::fast + solve_opts::no_approx);
  }

  res1 = Y - X * B1;
  while ((diff > conv) & (iter < maxiter)) {
    B0 = B1;
    for (j = 0; j < p; j++) {
      res1j = res1 +  X.col(j)* B1.row(j); //n q
      XRj =   trans(X.col(j)) * res1j;    //1 q
      rowvec t1=XRj/as_scalar(sh(j))*max(0.0,1-lam(j)/sqrt(accu(square(XRj))));
      B1.row(j) = t1;
      res1 = res1j - X.col(j)* B1.row(j);
    }
    l2B1 = accu(square(B1));
    if (l2B1 == 0) {
      iter = maxiter; 
    } else {
      diff = sqrt(accu(square(B0 - B1))/l2B1);
      iter = iter + 1;
    }
  }
  return(B1);
}



// [[Rcpp::export]]
Rcpp::List multicvrsolver(arma::mat Y, Rcpp::List Xlist, int rank, double eta, arma::vec Lam,
	std::string family, Rcpp::List Wini, std::string penalty, Rcpp::List opts){
	int K = Xlist.size();
	int n = Y.n_rows;
	int q = 1;
	int nrank = rank, maxIters = opts["maxIters"];
	bool  standardization = as<bool>(opts["standardization"]);
	double  tol = opts["tol"], h = 2, heps = 1;

	field <mat> X(K), XX(K), XW(K), W(K), B(K), Z1(K), Z2(K), Z3(K), Z4(K), GW(K), GZ1(K), GZ2(K), GZ3(K), GZ4(K), W0(K), C(K);

	arma::vec p(K), Wknorm, funcVal, diff = ones(1);
	arma::vec Wnz(K);   //Wnz[k]=1: W[k] is all 0
	uvec nzid;
	Rcpp::List obj;   //outputs

	arma::mat alpha1 = zeros(1, q), alpha2 = zeros(1, q), alpha3 = zeros(1, q), alpha4 = zeros(1, q), beta1 = zeros(nrank, q), beta2 = zeros(nrank, q), beta3 = zeros(nrank, q), beta4 = zeros(nrank, q);
	arma::mat sumCC = zeros(nrank + 1, nrank + 1);           // (r+1)-by-(r+1)
	arma::mat sumCZ1 = zeros(nrank + 1, q);                 // (r+1)-by-q
	arma::mat sumCZ2 = zeros(nrank + 1, q);
	arma::mat sumCZ3 = zeros(nrank + 1, q);
	arma::mat sumCZ4 = zeros(nrank + 1, q);
	arma::mat bhat1 = zeros(nrank + 1, q);
	arma::mat bhat2 = zeros(nrank + 1, q);
	arma::mat bhat3 = zeros(nrank + 1, q);
	arma::mat bhat4 = zeros(nrank + 1, q);
	arma::mat Y1 = zeros(n, 1);
	arma::mat Y2 = zeros(n, 1);
	arma::mat Y3 = zeros(n, 1);
	arma::mat Y4 = zeros(n, 1);
	arma::mat sumB = zeros(n, nrank);                   //n by r
	arma::mat Ztmp, resNorm;        
	arma::mat Mk, IrXk, vecW0k, BGWk, vecBGWk, onesXWk, IrXkalp, vecBGWk0;
	arma::mat  Gchol, U, V;
	arma:vec d, s;
	double  tmp1 = 0, tmp2 = 0, tmp3 = 0, t1, t2, difftmp;
 
	if (standardization){
		for (int k = 0; k < K; k++){
			X[k] = as<mat>(Xlist[k]) - ones(n)*mean(as<mat>(Xlist[k]));
			X[k] = X[k] / (ones(n)*stddev(as<mat>(Xlist[k])));
		}
	}
	else {
		for (int k = 0; k < K; k++){
			X[k] = as<mat>(Xlist[k]);
		}
	}

	for (int k = 0; k<K; k++){
		p[k] = X[k].n_cols;
		XX[k] = X[k].t()*X[k];
		W[k] = as<mat>(Wini[k]);
		XW[k] = X[k] * W[k];
		B[k] = XW[k];
		GW[k] = zeros(n, nrank);
		Z1[k] = zeros(n, q);
		Z2[k] = zeros(n, q);
		Z3[k] = zeros(n, q);
		Z4[k] = zeros(n, q);
		GZ1[k] = zeros(n, q);
		GZ2[k] = zeros(n, q);
		GZ3[k] = zeros(n, q);
		GZ4[k] = zeros(n, q);
	}
	//initialize Y1,Y2,Y3,Y4
	for (int i = 0; i < n; i++)
	{
		if (Y(i, 0) == 1)
			Y1(i, 0) = 1;
		else
			Y1(i, 0) = 0;
	}

	for (int i = 1; i < n; i++)
	{
		if (Y(i, 0) == 2)
			Y2(i, 0) = 1;
		else
			Y2(i, 0) = 0;
	}

	for (int i = 1; i < n; i++)
	{
		if (Y(i, 0) == 3)
			Y3(i, 0) = 1;
		else
			Y3(i, 0) = 0;
	}

	for (int i = 1; i < n; i++)
	{
		if (Y(i, 0) == 4)
			Y4(i, 0) = 1;
		else
			Y4(i, 0) = 0;
	}

	if (family == "multinomial")
	{
		for (int k = 0; k < K; k++)
		{
			Z1[k] = Newtonmultinomial(GZ1[k], B[k], alpha1, beta1, 1 - eta, Y1, h, Z1[k], Z2[k], Z3[k], Z4[k]);
			Z2[k] = Newtonmultinomial(GZ2[k], B[k], alpha2, beta2, 1 - eta, Y2, h, Z2[k], Z1[k], Z3[k], Z4[k]);
			Z3[k] = Newtonmultinomial(GZ3[k], B[k], alpha3, beta3, 1 - eta, Y3, h, Z3[k], Z1[k], Z2[k], Z4[k]);
			Z4[k] = Newtonmultinomial(GZ4[k], B[k], alpha4, beta4, 1 - eta, Y4, h, Z4[k], Z1[k], Z2[k], Z3[k]);
		}
	}

	//iteration begins
	int iter = 0;

	for (iter = 0; iter < maxIters; iter++){
		Wnz = zeros(K);
		sumCC = zeros(nrank + 1, nrank + 1);
		sumCZ1 = zeros(nrank + 1, q);
		sumCZ2 = zeros(nrank + 1, q);
		sumCZ3 = zeros(nrank + 1, q);
		sumCZ4 = zeros(nrank + 1, q);
		sumB = zeros(n, nrank);
		for (int k = 0; k < K; k++){
			W0[k] = W[k];
			if (accu(abs(W[k])) == 0){ Wnz[k] = 1; }
		}

		if (accu(Wnz) != 0) break;      // some Wk is all 0

		// beta-step
		for (int k = 0; k < K; k++){
			C[k] = join_rows(ones(n, 1), B[k]);                    // n by r+1
			sumCC = sumCC + C[k].t()*C[k];
			sumCZ1 = sumCZ1 + C[k].t()*(Z1[k] + GZ1[k] / h);
		}

		bhat1 = solve(sumCC + diagmat(0.00000001 * ones(nrank + 1)), sumCZ1);
		alpha1 = bhat1.row(0);                         // 1 by q    
		beta1 = bhat1.rows(1, nrank);                  // r by q


		for (int k = 0; k < K; k++){
			sumCZ2 = sumCZ2 + C[k].t()*(Z2[k] + GZ2[k] / h);
		}

		bhat2 = solve(sumCC + diagmat(0.00000001 * ones(nrank + 1)), sumCZ2);
		alpha2 = bhat2.row(0);                         // 1 by q    
		beta2 = bhat2.rows(1, nrank);                  // r by q


		for (int k = 0; k < K; k++){
			sumCZ3 = sumCZ3 + C[k].t()*(Z3[k] + GZ3[k] / h);
		}

		bhat3 = solve(sumCC + diagmat(0.00000001 * ones(nrank + 1)), sumCZ3);
		alpha3 = bhat3.row(0);                         // 1 by q    
		beta3 = bhat3.rows(1, nrank);                  // r by q

		for (int k = 0; k < K; k++){
			sumCZ4 = sumCZ4 + C[k].t()*(Z4[k] + GZ4[k] / h);
		}

		bhat4 = solve(sumCC + diagmat(0.00000001 * ones(nrank + 1)), sumCZ4);
		alpha4 = bhat4.row(0);                         // 1 by q    
		beta4 = bhat4.rows(1, nrank);                  // r by q


		// B-step
		for (int k = 0; k < K; k++){
			sumB = zeros(n, nrank);
			for (int j = 0; j < K; j++){
				sumB = sumB + B[j];
			}
			Mk = eta / 2 * (sumB - B[k]) + h / 2 * (XW[k] + GW[k] / h) + h / 2 * (Z1[k] - ones(n, 1)*alpha1 + GZ1[k] / h)*beta1.t() + h / 2 * (Z2[k] - ones(n, 1)*alpha2 + GZ2[k] / h)*beta2.t() + h / 2 * (Z3[k] - ones(n, 1)*alpha3 + GZ3[k] / h)*beta3.t() + h / 2 * (Z4[k] - ones(n, 1)*alpha4 + GZ4[k] / h)*beta4.t();

			bool success = svd_econ(U, s, V, Mk);
			if (!success) {
				Rcpp::Rcout << "OrthProc failed!" << std::endl;
			}
			B[k] = U*V.t();

		}

		// W-step
		for (int k = 0; k < K; k++){
			BGWk = B[k] - GW[k] / h;

			if (penalty == "GL1"){    //row-wise lasso penalty 
				Wknorm = sum(square(W0[k]), 1);
				//remove zero columns of X
				nzid = find(Wknorm != 0);
				W[k].rows(nzid) = MGlasso_C(BGWk, X[k].cols(nzid), ones(nzid.size())*Lam[k] / h, W0[k].rows(nzid), 1e-2, 30);
				//MGltmp = MGlasso_C(BGWk, X[k].cols(nzid), ones(nzid.size())*Lam[k]/h, W0[k].rows(nzid), 1e-2, 50);
				//W[k].rows(nzid)= as<mat>(MGltmp["B"]);
			}
			else if (penalty == "L1"){
				//entrywise lasso penalty: vectorize, then lasso
				vecW0k = reshape(W0[k], p[k] * nrank, 1);
				IrXk = kron(eye(nrank, nrank), X[k]);
				vecBGWk = reshape(BGWk, n*nrank, 1);
				//MGltmp = MGlasso_C(vecBGWk, IrXk, ones(p[k]*nrank)*Lam[k]/h, vecW0k, 1e-2, 50);
				//vecW0k = as<mat>(MGltmp["B"]);
				vecW0k = MGlasso_C(vecBGWk, IrXk, ones(p[k] * nrank)*Lam[k] / h, vecW0k, 1e-2, 30);
				W[k] = reshape(vecW0k, p[k], nrank);
			}
			else if (penalty == "enet"){
				//entrywise enet penalty: vectorize, then enet, alp=0.1(L2 penalty)
				vecW0k = reshape(W0[k], p[k] * nrank, 1);
				IrXk = kron(eye(nrank, nrank), X[k]);
				double alp = 0.1;
				IrXkalp = join_cols(IrXk, sqrt(alp)*eye(p[k] * nrank, p[k] * nrank));
				vecBGWk = reshape(BGWk, n*nrank, 1);
				vecBGWk0 = join_cols(vecBGWk, zeros(p[k] * nrank, 1));
				//MGltmp = MGlasso_C(vecBGWk0, IrXkalp, ones(p[k]*nrank)*Lam[k]*(1-alp)/h, vecW0k, 1e-2, 50);
				//vecW0k = as<mat>(MGltmp["B"]);
				vecW0k = MGlasso_C(vecBGWk0, IrXkalp, ones(p[k] * nrank)*Lam[k] * (1 - alp) / h, vecW0k, 1e-2, 30);
				vecW0k = (1 + alp*Lam[k])*vecW0k;
				W[k] = reshape(vecW0k, p[k], nrank);
			}

			XW[k] = X[k] * W[k];    // store XW[k] fo later use
		}

		// Z-step 
		for (int k = 0; k < K; k++)
		{
			Z1[k] = Newtonmultinomial(GZ1[k], B[k], alpha1, beta1, 1 - eta, Y1, h, Z1[k], Z2[k], Z3[k], Z4[k]);
			Z2[k] = Newtonmultinomial(GZ2[k], B[k], alpha2, beta2, 1 - eta, Y2, h, Z2[k], Z1[k], Z3[k], Z4[k]);
			Z3[k] = Newtonmultinomial(GZ3[k], B[k], alpha3, beta3, 1 - eta, Y3, h, Z3[k], Z1[k], Z2[k], Z4[k]);
			Z4[k] = Newtonmultinomial(GZ4[k], B[k], alpha4, beta4, 1 - eta, Y4, h, Z4[k], Z1[k], Z2[k], Z3[k]);
		}

		// Dual-step 
		for (int k = 0; k < K; k++){
			GW[k] = GW[k] + h*(XW[k] - B[k]);
			GZ1[k] = GZ1[k] + h*(Z1[k] - B[k] * beta1 - ones(n, 1)*alpha1);
			GZ2[k] = GZ2[k] + h*(Z2[k] - B[k] * beta2 - ones(n, 1)*alpha2);
			GZ3[k] = GZ3[k] + h*(Z3[k] - B[k] * beta3 - ones(n, 1)*alpha3);
			GZ4[k] = GZ4[k] + h*(Z4[k] - B[k] * beta4 - ones(n, 1)*alpha4);
		}

		// objective function
		tmp1 = 0;
		for (int k = 1; k < K; k++){
			for (int j = 0; j < k; j++){
				tmp1 = tmp1 + accu(square(XW[k] - XW[j])) / 2;
			}
		}

		tmp2 = 0;
		for (int k = 0; k < K; k++)
		{
			tmp2 = tmp2 - accu(Y1%Z1[k] + Y2%Z2[k] + Y3%Z3[k] + Y4%Z4[k] - log(exp(Z1[k] + Z2[k] + Z3[k] + Z4[k]))) / n;
		}

		tmp3 = 0;
		for (int k = 0; k < K; k++){
			if (penalty == "GL1"){
				tmp3 = tmp3 + Lam[k] * accu(sqrt(sum(square(W[k]), 1)));
			}
			else if (penalty == "L1"){
				tmp3 = tmp3 + Lam[k] * accu(abs(W[k]));
			}
		}

		funcVal = join_cols(funcVal, (eta*tmp1 + (1 - eta)*tmp2 + tmp3)*ones(1, 1));

		if (iter > 0){
			t1 = as_scalar(funcVal(iter) - funcVal(iter - 1));
			t2 = abs(as_scalar(funcVal(iter - 1)));
			difftmp = t1 / t2;
			diff = join_cols(diff, difftmp*ones(1, 1));
			if (abs(difftmp) < tol) break;
		}

		h = h*heps;
		if (iter == maxIters - 1){
			Rcpp::Rcout << "Not converge yet!" << difftmp << std::endl;
		}
	}     //end of for
  
  
  
  //obj["family"] = family;
	obj["iter"] = iter;
	//obj["X"] = X;
	obj["W"] = W;
	obj["B"] = B;
	obj["Z1"] = Z1;
	obj["Z2"] = Z2;
	obj["Z3"] = Z3;
	obj["Z4"] = Z4;
	obj["alpha1"] = alpha1;
	obj["alpha2"] = alpha2;
	obj["alpha3"] = alpha3;
	obj["alpha4"] = alpha4;
	obj["beta1"] = beta1;
	obj["beta2"] = beta2;
	obj["beta3"] = beta3;
	obj["beta4"] = beta4;
	obj["objvals"] = funcVal;
	//obj["diff"] = diff;
	obj["ver"] = "test";
	return(obj);
}


