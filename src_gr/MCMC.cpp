#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class Session {
public:
  int numSamples;
  int numKernels; // TODO test it's only 1/2
  double maternCoef;
  double robbinsDelta;
  double robbinsB;
  
  Session(Rcpp::List sessionSettings){
    numSamples = Rcpp::as<int>(sessionSettings["numSamples"]);
    numKernels = Rcpp::as<int>(sessionSettings["numKernels"]);
    maternCoef = Rcpp::as<double>(sessionSettings["maternCoef"]);
    robbinsDelta = Rcpp::as<double>(sessionSettings["robbinsDelta"]);
    robbinsB = Rcpp::as<double>(sessionSettings["robbinsB"]);
  }
  
};

class Data {
public:
  int numDays;
  int numGroups;
  Rcpp::NumericVector t; // should be ordered and equally spaced
  Rcpp::NumericVector pos;
  Rcpp::NumericVector tot;
  Rcpp::NumericVector dayToGroup;
  double mEta;
  double sigEta;
  double aW;
  double bW;
  Rcpp::NumericVector ker_logTheta0; // { logSigma_K1, logRange_K1, logSigma_K2, logRange_K2 }
  arma::mat ker_B;
  arma::mat ker_Binv;
  
  Data(Rcpp::List inputData){
    numDays = Rcpp::as<int>(inputData["num_days"]);
    numGroups = Rcpp::as<int>(inputData["num_groups"]);
    t = Rcpp::as<Rcpp::NumericVector>(inputData["t"]);
    pos = Rcpp::as<Rcpp::NumericVector>(inputData["pos"]);
    tot = Rcpp::as<Rcpp::NumericVector>(inputData["tot"]);
    dayToGroup = Rcpp::as<Rcpp::NumericVector>(inputData["day_to_group"]);
    mEta = Rcpp::as<double>(inputData["m_eta"]);
    sigEta = Rcpp::as<double>(inputData["sig_eta"]);
    aW = Rcpp::as<double>(inputData["a_w"]);
    bW = Rcpp::as<double>(inputData["b_w"]);
    ker_logTheta0 = Rcpp::as<Rcpp::NumericVector>(inputData["log_theta_0"]); // { logSigma_K1, logRange_K1, logSigma_K2, logRange_K2 }
    ker_B = Rcpp::as<arma::mat>(inputData["B"]); // nxn matrix with n = size of ker_logTheta0 // TODO test
    ker_Binv = ker_B.i();
  }
  
};

class Covariance {
public:
  arma::vec covarianceMatrix;
  //arma::mat inverseCovariance;
  double *choleskyVecCovariance;
  arma::mat choleskyMatCovariance;
  arma::vec u;
  double uTu;
  double logDetCovariance;
  //double detCovariance;
  
  Covariance() {
    covarianceMatrix = arma::vec();
    choleskyMatCovariance = arma::mat();
    u = arma::vec();
  }
  
  Covariance(int size) {
    covarianceMatrix = arma::vec(size);
    choleskyVecCovariance = new double[size*size];
    choleskyMatCovariance = arma::mat(size, size);
    u = arma::vec(size);
  }
  
  void computeCovarianceMatrix(arma::vec &distanceMatrix, double maternCoef, int numKernels,
                               arma::vec &covHyperparam) {
    // TODO check maternCoef is 0.5 or 1.5
    // TODO compute only the first row and use arma::toeplitz
    
    double ker1_sigma = exp(covHyperparam(0));
    double ker1_length = exp(covHyperparam(1))/2;
    double ker2_sigma = 0;
    //double ker2_length = 0;
    double ker2_period = 0;
    if (numKernels == 2) {
      ker2_sigma = exp(covHyperparam(2));
      //ker2_length = exp(covHyperparam(3))/2;
      ker2_period = exp(covHyperparam(3))/2;
    }
    
    // Loop
    for (int i = 0; i < covarianceMatrix.n_rows; i++) {
      // Matern kernel
      if (maternCoef == 0.5) {
        covarianceMatrix(i) = ker1_sigma*exp(-distanceMatrix(i)/ker1_length);
      } else if (maternCoef == 1.5) {
        covarianceMatrix(i) = ker1_sigma*
          (1 + (sqrt(3)*distanceMatrix(i)/ker1_length))*exp(-sqrt(3)*distanceMatrix(i)/ker1_length);
      }
      // 
      // Periodic kernel
      if (numKernels == 2) {
        //covarianceMatrix(i) += ker2_sigma*exp(-2*sin(distanceMatrix(i)/2)*sin(distanceMatrix(i)/2)/(ker2_length*ker2_length));
        covarianceMatrix(i) += ker2_sigma*exp( -distanceMatrix(i)*distanceMatrix(i)/(2*3600)
                                                 -2*sin(3.141593*distanceMatrix(i)/ker2_period)*sin(3.141593*distanceMatrix(i)/ker2_period)  );
      }
    }
  }
  
  //void computeInverseMatrix(arma::mat &inverseCovariance, arma::mat covarianceMatrix) {
  //  inverseCovariance = arma::inv_sympd(covarianceMatrix);
  //}
  
  //void computeCholeskyArmaMatrix() {
  //  // compute LL decomposition using armadillo
  //  choleskyCovariance = arma::chol(covarianceMatrix, "lower");
  //}
  
  void computeCholeskyToeplitzMatrix() {
    // compute LL decomposition using algorithm for Toeplitz matrices
    
    // l -> choleskyTCovariance, t -> covarianceMatrix <- t, n -> covarianceMatrix.n_rows
    int n = covarianceMatrix.n_rows;
    double div;
    double *g;
    double g1j;
    double g2j;
    int i;
    int j;
    double rho;
    
    g = ( double * ) malloc ( 2 * n * sizeof ( double ) );
    
    for ( j = 0; j < n; j++ )
    {
      g[0+j*2] = covarianceMatrix[j]; //t[j];
    }
    g[1+0*2] = 0.0;
    for ( j = 1; j < n; j++ )
    {
      g[1+j*2] = covarianceMatrix[j]; //t[j];
    }
    
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        choleskyVecCovariance[i+j*n] = 0.0;
      }
    }
    
    for ( i = 0; i < n; i++ )
    {
      choleskyVecCovariance[i+0*n] = g[0+i*2];
    }
    
    for ( j = n - 1; 1 <= j; j-- )
    {
      g[0+j*2] = g[0+(j-1)*2];
    }
    
    g[0+0*2] = 0.0;
    
    for ( i = 1; i < n; i++ )
    {
      rho = - g[1+i*2] / g[0+i*2];
      div = sqrt ( ( 1.0 - rho ) * ( 1.0 + rho ) );
      for ( j = i; j < n; j++ )
      {
        g1j = g[0+j*2];
        g2j = g[1+j*2];
        g[0+j*2] = (       g1j + rho * g2j ) / div;
        g[1+j*2] = ( rho * g1j +       g2j ) / div;
      }
      for ( j = i; j < n; j++ )
      {
        choleskyVecCovariance[j+i*n] = g[0+j*2];
      }
      for ( j = n - 1; i < j; j-- )
      {
        g[0+j*2] = g[0+(j-1)*2];
      }
      g[0+i*2] = 0.0;
    }
    
    delete [] g;
    
    // choleskyMatCovariance = arma::toeplitz(choleskyVecCovariance); // doesnt work for pointers
    
    // TODO is this efficient?
    // we need it only for using 'solve' in computeXSX()
    // TODO maybe easier solving it ourselves?
    for (int i = 0; i < covarianceMatrix.n_rows; i++) {
      for (int j = 0; j <= i; j++) {
        choleskyMatCovariance(i,j) = choleskyVecCovariance[i + (j*covarianceMatrix.n_rows)];//filling by cols!
      }
    }
  }
  
  void computeDetFromCholesky(){
    // compute determinant from L
    logDetCovariance = 0;
    for (int i = 0; i < covarianceMatrix.n_rows; i++) {
      //logDetTimes2 += log(choleskyCovariance(i,i));
      logDetCovariance += 2*log(choleskyVecCovariance[i + i*covarianceMatrix.n_rows]);
    }
  }
  
  void computeXSX(arma::vec &x) {
    // compute xT*inverse*x = z*z
    u = arma::solve(arma::trimatl(choleskyMatCovariance), x);
    uTu = arma::as_scalar(u.t()*u);
  }
};

class State {
public:
  double logEta;
  double tauW;
  arma::vec w;
  arma::vec x;
  arma::mat sigmaWMinus;//TODO should it be somewhere else?
  arma::vec distanceMatrix;//TODO should it be somewhere else?
  arma::vec covHyperparam;  // { logSigma_K1, logRange_K1, logSigma_K2, logRange_K2 } // TODO should I move it to Covariance?
  Covariance covarianceClass;
  
  State(Data &data, Session &session) {
    w = arma::vec(data.numGroups);
    x = arma::vec(data.numDays);
    distanceMatrix = arma::vec(data.numDays);
    
    // Construct sigmaWMinus matrix. See notes 21.05.2024
    sigmaWMinus = arma::mat(data.numGroups, data.numGroups);
    sigmaWMinus(0,0) = 6.0;
    for(int i = 1; i < data.numGroups; i++){
      sigmaWMinus(i, 0) = -1.0;
      sigmaWMinus(0, i) = -1.0;
      sigmaWMinus(i, i) = 1.0;
    }
    
    // Construct distanceMatrix
    // Since equally spaced, we store less
    for (int i = 0; i < data.numDays; i++) {
      distanceMatrix(i) = abs(data.t[i] - data.t[0]);
    }
    
    covHyperparam = arma::vec(2*session.numKernels);
    covarianceClass = Covariance(data.numDays);
  }
  
  void initialiseStateRandom(Data &data, Session &session){
    logEta = R::rnorm(data.mEta, data.sigEta);
    tauW = R::rgamma(data.aW, 1/data.bW);
    //ker1_logRange = R::rnorm(data.ker1_logTheta0[0], data.ker1_B(0,0));
    //ker1_logSigma = R::rnorm(data.ker1_logTheta0[1] + data.ker1_B(0,1)*(ker1_logRange - data.ker1_logTheta0[0])/data.ker1_B(0,0),
    //                    pow((data.ker1_B(1,1) - pow(data.ker1_B(0,1),2))/data.ker1_B(0,0), 0.5));
    //ker2_logRange = R::rnorm(data.ker2_logTheta0[0], data.ker2_B(0,0));
    //ker2_logSigma = R::rnorm(data.ker2_logTheta0[1] + data.ker2_B(0,1)*(ker2_logRange - data.ker2_logTheta0[0])/data.ker2_B(0,0),
    //                         pow((data.ker2_B(1,1) - pow(data.ker2_B(0,1),2))/data.ker2_B(0,0), 0.5));
    w = Rcpp::rnorm(w.n_elem, 0, pow(data.bW/data.aW, 0.5));
    x = Rcpp::rnorm(x.n_elem, 0, exp(data.ker_logTheta0[1]));
    covHyperparam = arma::vec(2*session.numKernels);
    covHyperparam(0) = R::rnorm(data.ker_logTheta0[0], data.ker_B(0,0));
    covHyperparam(1) = R::rnorm(data.ker_logTheta0[1], data.ker_B(1,1));
    if(session.numKernels == 2){
      covHyperparam(2) = R::rnorm(data.ker_logTheta0[2], data.ker_B(2,2));
      covHyperparam(3) = R::rnorm(data.ker_logTheta0[3], data.ker_B(3,3));
    }
  }
  
  void initialiseFromRList(Rcpp::List &list, Session &session) {
    logEta = Rcpp::as<double>(list["logEta"]);
    tauW = Rcpp::as<double>(list["tauW"]);
    w = Rcpp::as<Rcpp::NumericVector>(list["w"]);
    x = Rcpp::as<Rcpp::NumericVector>(list["x"]);
    covHyperparam = Rcpp::as<Rcpp::NumericVector>(list["log_theta"]); // { logSigma_K1, logRange_K1, logSigma_K2, logRange_K2 }
  }
  
};

class Cache {
public:
  double logEta;
  Rcpp::NumericVector w;
  Rcpp::NumericVector x;
  Rcpp::NumericVector covHyperparam;
  
  Cache(int numDays, int numGroups, int numKernels, double value) {
    logEta = value;
    w = Rcpp::NumericVector(numGroups, value);
    x = Rcpp::NumericVector(numDays, value);
    covHyperparam = Rcpp::NumericVector(2*numKernels, value);
  }
};

class Storage {
public:
  arma::vec logEta;
  arma::vec tauW;
  arma::mat w;
  arma::mat x;
  arma::mat covHyperparam;
  
  Storage(int samples, int numDays, int numGroups, int numKernels) {
    logEta = arma::vec(samples);
    tauW = Rcpp::NumericVector(samples);
    w = arma::mat(numGroups, samples);
    x = arma::mat(numDays, samples);
    covHyperparam = arma::mat(2*numKernels, samples);
  }
};

double getLogDensityEta(double logEta, Data &data, State &state) {
  double eta = exp(logEta);
  double prior = -0.5*(logEta - data.mEta)*(logEta - data.mEta)/(data.sigEta*data.sigEta);
  double likelihood = 0;
  for (int i = 0; i < data.numDays; i++) {
    double muT = exp(state.x(i) + state.w(data.dayToGroup(i) - 1));
    likelihood += data.pos[i]*log(muT/eta) - (eta + data.pos[i])*log(1 + (muT/eta));
  }
  return prior + likelihood;
}

void updateEta(double it, Data &data, Session &session, State &state, Cache &jumps, Cache &acceptance, Cache &rejection) {
  double proposal = R::rnorm(state.logEta, jumps.logEta);
  double lar = getLogDensityEta(proposal, data, state) - getLogDensityEta(state.logEta, data, state);
  double u = R::runif(0,1);
  if(log(u) < lar) {
    state.logEta = proposal;
    acceptance.logEta += 1;
    //std::cout << "accepted ";
  }else{
    rejection.logEta += 1;
    //std::cout << "rejected ";
  }
  jumps.logEta *= exp(session.robbinsDelta*(std::min<double>(1, exp(lar)) - 0.44)/std::max<double>(it - session.robbinsB, 1));
}

void updateTauW(Data &data, State &state) {
  double wSw = 0;
  for (int i = 0; i < state.sigmaWMinus.n_rows; i++) {
    for (int j = 0; j < state.sigmaWMinus.n_rows; j++) {
      wSw += state.w(i)*state.sigmaWMinus(i,j)*state.w(j);
    }
  }
  state.tauW = R::rgamma(data.aW + (data.numGroups - 1)/2, 1/(data.bW + 0.5*wSw));
}

double getLogDensityW(arma::vec &w, Data &data, State &state) {
  // See notes 21.05.2024
  double prior = 0;
  for (int j = 0; j < data.numGroups; j++) {
    for (int k = 0; k < data.numGroups; k++) {
    prior += w(j)*state.sigmaWMinus(j,k)*w(k);
    }
  }
  prior *= -0.5*state.tauW;
  double likelihood = 0;
  for (int i = 0; i < data.numDays; i++) {
    double muT = exp(state.x(i) + w(data.dayToGroup(i) - 1));
    likelihood += data.pos[i]*log(muT/exp(state.logEta)) - (exp(state.logEta) + data.pos[i])*log(1 + (muT/exp(state.logEta)));
  }
  return prior + likelihood;
}

void updateW(double it, Data &data, Session &session, State &state, Cache &jumps, Cache &acceptance, Cache &rejection) {
  // See notes 21.05.2024 (end bit)
  // Single site updates, adjusting all before computing log-density
  arma::vec proposalW(state.w.n_elem);
  for (int i = 0; i < state.w.n_elem; i++) {
    proposalW = state.w;
    proposalW(i) = R::rnorm(state.w(i), jumps.w(i));
    proposalW = proposalW - mean(proposalW);
    double lar = getLogDensityW(proposalW, data, state) - getLogDensityW(state.w, data, state);
    double u = R::runif(0,1);
    if(log(u) < lar) {
      state.w = proposalW;
      acceptance.w(i) += 1;
      //if(i == 0) {
      //  std::cout << "accepted ";
      //}
    }else{
      rejection.w(i) += 1;
      //if(i == 0) {
      //  std::cout << "rejected ";
      //}
    }
    jumps.w(i) *= exp(session.robbinsDelta*(std::min<double>(1, exp(lar)) - 0.44)/std::max<double>(it - session.robbinsB, 1));
  }
}

double getLogDensityXi(int i, double xi, double uTu, Data &data, State &state) {
  // See notes 21.05.2024
  //double prior = -0.5*uTu;
  //double muT = exp(xi + state.w(data.dayToGroup(i) - 1));
  //double likelihood = data.pos[i]*log(muT/exp(state.logEta)) - (exp(state.logEta) + data.pos[i])*log(1 + (muT/exp(state.logEta)));
  //return prior + likelihood;
  double muT_eta = exp(xi + state.w(data.dayToGroup(i) - 1))/exp(state.logEta);
  double priorLikelihood = -0.5*uTu + data.pos[i]*log(muT_eta) - (exp(state.logEta) + data.pos[i])*log(1 + muT_eta);
  return priorLikelihood;
}

void updateX(double it, Data &data, Session &session, State &state, Cache &jumps, Cache &acceptance, Cache &rejection) {//arma::mat &auxLuXupdate
  // See notes 22.05.2024 and picture of board
  // Single site updates, adjusting all before computing log-density
  
  // Construct auxLuXupdate
  // auxLuXupdate is Lu but with a cumSum from left to right
  //for (int i = 0; i < data.numDays; i++) {
  //  auxLuXupdate(i,0) = state.covarianceClass.choleskyVecCovariance[i + (0*data.numDays)]*state.covarianceClass.u[0];
  //}
  //for (int j = 1; j < data.numDays - 1; j++) {
  //  for (int i = j; i < data.numDays; i++) {
  //    auxLuXupdate(i,j) = auxLuXupdate(i,j - 1) + state.covarianceClass.choleskyVecCovariance[i + (j*data.numDays)]*state.covarianceClass.u[j];
  //  }
  //}
  
  // Update x - single update, from n to 1
  arma::vec proposalX = state.x;
  arma::vec proposal_u;// = state.covarianceClass.u;
  double proposal_uTu;
  double lar;
  double u;
  //arma::vec v;
  for (int j = data.numDays - 1; j >= 0; j--) {
    // Update xi
    proposalX(j) = R::rnorm(state.x(j), jumps.x(j));
    
    // Compute u
    // Opt. 1 - my algorithm using auxLuXupdate
    //if (j > 0) {
    //  v = arma::solve(arma::trimatl(state.covarianceClass.choleskyMatCovariance.submat(j, j, data.numDays - 1, data.numDays - 1)),
    //                  state.x.subvec(j, data.numDays - 1) - auxLuXupdate.submat(j, j - 1, data.numDays - 1, j - 1), arma::solve_opts::fast);
    //  proposal_u.subvec(j, data.numDays - 1) = v;
    //} else {
    //  proposal_u = arma::solve(arma::trimatl(state.covarianceClass.choleskyMatCovariance), proposalX, arma::solve_opts::fast);
    //}
    // Opt. 2 - solving the system each time (slightly more efficient)... 25 min for 5000 it
    //proposal_u = arma::solve(arma::trimatl(state.covarianceClass.choleskyMatCovariance), proposalX); // inefficient
    proposal_u = arma::solve(arma::trimatl(state.covarianceClass.choleskyMatCovariance), proposalX, arma::solve_opts::fast); // better
    //proposal_u = state.covarianceClass.u; // TODO BORRAR. Only for test speed
    
    // Compute uTu
    proposal_uTu = arma::as_scalar(proposal_u.t()*proposal_u);// = 16 sec for 5000 iterations
    //proposal_uTu = state.covarianceClass.uTu; // TODO BORRAR. Only for test speed
    
    // Accept/reject
    lar = getLogDensityXi(j, proposalX(j), proposal_uTu, data, state) - getLogDensityXi(j, state.x(j), state.covarianceClass.uTu, data, state);
    u = R::runif(0,1);
    if(log(u) < lar) {
      state.x(j) = proposalX(j);
      acceptance.x(j) += 1;
      state.covarianceClass.u = proposal_u;
      state.covarianceClass.uTu = proposal_uTu;
    }else{
      rejection.x(j) += 1;
    }
    jumps.x(j) *= exp(session.robbinsDelta*(std::min<double>(1, exp(lar)) - 0.44)/std::max<double>(it - session.robbinsB, 1));
  }
}

double getLogDensityCovHyperparam(arma::vec &hyperparameters, Covariance &covariance, Data &data, State &state) {
  // Block update of all hyperparameters
  // hyperparameters: { K1 sigma , K1 range , K2 sigma , K2 period }
  double prior = -0.5*arma::as_scalar(hyperparameters.t()*data.ker_Binv*hyperparameters);
  double likelihood = -0.5*covariance.logDetCovariance - 0.5*covariance.uTu;
  
  return prior + likelihood;
}

void updateCovHyperparam(int it, Covariance &proposalClass, Data &data, Session &session, State &state, Cache &jumps, Cache &acceptance, Cache &rejection) {
  // Block update of all hyperparameters. Acceptance/rejection/jumps are stored on the first entrance
  arma::vec proposal(4); // { K1 sigma , K1 range , K2 sigma , K2 range }
  
  // TODO propose block update with Munroe adaptive
  // Individual jumps
  proposal(0) = R::rnorm(state.covHyperparam(0), jumps.covHyperparam(0));
  proposal(1) = R::rnorm(state.covHyperparam(1), jumps.covHyperparam(1));
  if(session.numKernels == 2) {
    proposal(2) = R::rnorm(state.covHyperparam(2), jumps.covHyperparam(2));
    proposal(3) = R::rnorm(state.covHyperparam(3), jumps.covHyperparam(3));
  }
  
  proposalClass.computeCovarianceMatrix(state.distanceMatrix, session.maternCoef, session.numKernels, proposal);
  proposalClass.computeCholeskyToeplitzMatrix();
  proposalClass.computeDetFromCholesky();
  proposalClass.computeXSX(state.x);
    
  double lar = getLogDensityCovHyperparam(proposal, proposalClass, data, state) -
    getLogDensityCovHyperparam(state.covHyperparam, state.covarianceClass, data, state);
    
    double u = R::runif(0,1);
    if(log(u) < lar) {
      state.covHyperparam(0) = proposal(0);
      state.covHyperparam(1) = proposal(1);
      if(session.numKernels == 2) {
        state.covHyperparam(2) = proposal(2);
        state.covHyperparam(3) = proposal(3);
      }
      // Stored on the first entry
      acceptance.covHyperparam(0) += 1;
    }else{
      rejection.covHyperparam(0) += 1;
    }
    jumps.covHyperparam(0) *= exp(session.robbinsDelta*(std::min<double>(1, exp(lar)) - 0.44)/std::max<double>(it - session.robbinsB, 1));
}

Rcpp::List createOutput(Storage &storage, Cache &jumps, Cache &acceptance, Cache &rejection, Session &session){
  Rcpp::List storageR(5);
  storageR(0) = storage.logEta;
  storageR(1) = storage.tauW;
  storageR(2) = storage.w;
  storageR(3) = storage.x;
  storageR(4) = storage.covHyperparam;
                                    
  Rcpp::List jumpsR(4);
  jumpsR(0) = jumps.logEta;
  jumpsR(1) = jumps.w;
  jumpsR(2) = jumps.x;
  jumpsR(3) = jumps.covHyperparam;
  
  Rcpp::List acceptanceR(4);
  acceptanceR(0) = acceptance.logEta;
  acceptanceR(1) = acceptance.w;
  acceptanceR(2) = acceptance.x;
  acceptanceR(3) = acceptance.covHyperparam;
  
  Rcpp::List rejectionR(4);
  rejectionR(0) = rejection.logEta;
  rejectionR(1) = rejection.w;
  rejectionR(2) = rejection.x;
  rejectionR(3) = rejection.covHyperparam;
  
  Rcpp::List sessionList(5);
  sessionList(0) = session.numSamples;
  sessionList(1) = session.numKernels;
  sessionList(2) = session.maternCoef;
  sessionList(3) = session.robbinsDelta;
  sessionList(4) = session.robbinsB;
  
  std::vector<std::string> namesStorage = {"logEta", "tauW", "w", "x", "covHyperparam"};
  std::vector<std::string> namesCache = {"logEta", "w", "x", "covHyperparam"};
  std::vector<std::string> namesSession = {"numSamples", "numKernels", "maternCoef", "robbinsDelta", "robbinsB"};
  storageR.attr("names") = Rcpp::wrap(namesStorage);
  jumpsR.attr("names") = Rcpp::wrap(namesCache);
  acceptanceR.attr("names") = Rcpp::wrap(namesCache);
  rejectionR.attr("names") = Rcpp::wrap(namesCache);
  sessionList.attr("names") = Rcpp::wrap(namesSession);
  
  Rcpp::List outputR(5);
  outputR(0) = storageR;
  outputR(1) = jumpsR;
  outputR(2) = acceptanceR;
  outputR(3) = rejectionR;
  outputR(4) = sessionList;
  std::vector<std::string> names = {"storage", "jumps", "acceptance", "rejection", "session"};
  outputR.attr("names") = Rcpp::wrap(names);
  
  return outputR;
}

// [[Rcpp::export]]
Rcpp::List runMCMC(Rcpp::List dataInput,
                   Rcpp::List sessionSettings,
                   bool provideInitialState,
                   Rcpp::List initialState,
                   Rcpp::Function determinantFn)  {
  // Input : same as input of rumModel(), numSamples
  // Output : same as samples in output of rumModel()
  
  // Create objects
  Data data(dataInput);
  Session session(sessionSettings);
  Storage storage(session.numSamples, data.numDays, data.numGroups, session.numKernels);
  State state(data, session); // includes auxiliar matrices
  Cache acceptance(data.numDays, data.numGroups, session.numKernels, 0);
  Cache rejection(data.numDays, data.numGroups, session.numKernels, 0);
  Cache jumps(data.numDays, data.numGroups, session.numKernels, 1);
  
  // Initialise values
  if(provideInitialState == true) {
    state.initialiseFromRList(initialState, session);
  } else {
    state.initialiseStateRandom(data, session);
  }
  
  // Compute delta matrix, inverse and determinant of delta matrix
  // TODO
  state.covarianceClass.computeCovarianceMatrix(state.distanceMatrix, session.maternCoef, session.numKernels, state.covHyperparam);
  state.covarianceClass.computeCholeskyToeplitzMatrix();
  state.covarianceClass.computeDetFromCholesky();
  state.covarianceClass.computeXSX(state.x);
  
  // Copy covarianceClass for proposals
  Covariance proposalCovarianceClass(data.numDays);
  proposalCovarianceClass = state.covarianceClass;
  
  // Create storage for cumulated Lu from left to right. For updateX
  //arma::mat auxLuXupdate(data.numDays, data.numDays);
  
  // Iterations
  std::cout << "Updating ";
  for (int it = 0; it < session.numSamples; it++) {
    std::cout << it << "... ";
    
    // Update
    updateEta(it, data, session, state, jumps, acceptance, rejection);
    updateTauW(data, state);
    updateW(it, data, session, state, jumps, acceptance, rejection);
    updateX(it, data, session, state, jumps, acceptance, rejection);
    // thetaKer1: monroe in block?    proposalCovarianceClass
    // thetaKer2: monroe in block?
    
    //state.covarianceClass.computeCovarianceMatrix(state.distanceMatrix, session.maternCoef, session.numKernels,
    //                                              state.ker1_logRange, state.ker1_logSigma, state.ker2_logRange, state.ker2_logSigma);
    //state.covarianceClass.computeCholeskyToeplitzMatrix();
    //state.covarianceClass.computeDetFromCholesky();
    //state.covarianceClass.computeXSX(state.x);
    
    // Store
    storage.logEta(it) = state.logEta;
    storage.tauW(it) = state.tauW;
    storage.w.col(it) = state.w;
    storage.x.col(it) = state.x;
    storage.covHyperparam.col(it) = state.covHyperparam;
  }
  
  std::cout << "Finished.\n";
  return createOutput(storage, jumps, acceptance, rejection, session);
}



