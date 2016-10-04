double muOrderZero(double ks){
    // limits were established with signa = 2.*3./sqrt(2.), limits should scale propto velDispersion/(2.*3./sqrt(2.))
    if(ks > 0.0914289)  return 0.5*pow(pi, 0.5)*gsl_sf_erf(ks)/ks;

    else{
        return 1.0 - 0.333333*pow(ks, 2.) + 0.1*pow(ks, 4.) - 0.0238095*pow(ks, 6.);
    }
}


double muOrderTwo(double ks){
    if(ks > 0.15727469)  return pow(4.*ks*ks*ks, -1.)*(pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-1.*ks*ks));
    
    else{
        return 1./3. - ks*ks/5. + pow(ks, 4.)/14.;
    }
}


double muOrderFour(double ks){
    if(ks>0.2418305)  return pow(8.*pow(ks, 5.), -1.)*(3.*pow(pi, 0.5)*gsl_sf_erf(ks) -6.*ks*exp(-1.*ks*ks) - 4.*pow(ks, 3.)*exp(-ks*ks));

    else{
        return 1./5. - ks*ks/7. + pow(ks, 4.)/18.;
    }
}


double muOrderSix(double ks){
    if(ks>0.335168)  return pow(16.*pow(ks, 7.), -1.)*(15.*pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-ks*ks)*(15. + 10.*ks*ks + 4.*pow(ks, 4.)));

    else{
         return 1./7. - ks*ks/9. + pow(ks, 4.)/22.;
    }
}

//Mike
double kaiserGauss_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}

//Mike
double kaiserGauss_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}
