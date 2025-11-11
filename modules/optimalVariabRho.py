import math 
import ROOT

def wVariab(genTauP4,genRhoP4,genPionP4,beamE, testAtau=0):

    # GEN 
    boostGEN=ROOT.TVector3()
    boostGEN=genTauP4.BoostVector()

    genRho_TauRes=ROOT.TLorentzVector()
    genRho_TauRes.SetPxPyPzE(genRhoP4.Px(),genRhoP4.Py(),genRhoP4.Pz(),genRhoP4.E())
    genRho_TauRes.Boost(-boostGEN)

    genBoostRho=ROOT.TVector3()
    genBoostRho=genRhoP4.BoostVector()

    genPion_RhoRes=ROOT.TLorentzVector()
    genPion_RhoRes.SetPxPyPzE(genPionP4.Px(),genPionP4.Py(),genPionP4.Pz(),genPionP4.E())
    genPion_RhoRes.Boost(-genBoostRho)

    v1 = genRho_TauRes.Vect()
    v2 = genTauP4.Vect()
    gen_theta_Rho= v1.Angle(v2)
    z = math.cos( gen_theta_Rho)
    #print (theta_Rho,z)

    x=genRhoP4.E()/genTauP4.E()
    sqrts=beamE*2

    mtau=1.7769
    mRho=genRhoP4.M()

    cos_theta=math.cos(gen_theta_Rho) 
    
    #cos_theta = (2 *x * mtau*mtau - mtau*mtau - mRho*mRho)/(mtau*mtau-mRho*mRho)/math.sqrt(1-4*mtau*mtau/sqrts/sqrts)
    #print ("%.2f %.2f  | %.2f  %.2f | %.2f %.2f | %.2f" %(genTauP4.Theta(),genRho_TauRes.Theta(),mRho, x, z,cos_theta, z/cos_theta))
    # From this printout: great agreement if Rho is really a Rho, bad if the mass is far from the Rho (==if the event is not well assigned)
    #gen_theta_Rho=math.acos(cos_theta)
    #z=cos_theta

    v1 = genPion_RhoRes.Vect()
    v2 = genRhoP4.Vect()
    gen_beta= v1.Angle(v2)
    gen_cosBeta = math.cos( gen_beta)
    #print (theta_Rho,z)

    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(x*x-4*mRho*mRho/sqrts/sqrts))
    if cosPsi>1:
       print ('What happened?', cosPsi)
       cosPsi=1
    if cosPsi<-1:
       print ('What happened?', cosPsi)
       cosPsi=-1
    anglePsi=math.acos(cosPsi)

    w_a= (-2+mtau*mtau/mRho/mRho + 2*(1+mtau*mtau/mRho/mRho)*(3*cosPsi-1)/2*(3*gen_cosBeta*gen_cosBeta-1)/2)* math.cos(gen_theta_Rho)
    w_b = 3 * math.sqrt(mtau*mtau/mRho/mRho) * (3*gen_cosBeta*gen_cosBeta-1)/2 * math.sin(2 * anglePsi) * math.sin(gen_theta_Rho)
    w_c = 2 + (mtau*mtau/mRho/mRho) - 2 *(1- (mtau*mtau/mRho/mRho)) * (3*cosPsi-1)/2 * (3*gen_cosBeta*gen_cosBeta-1)/2


    w= (w_a+w_b)/w_c

    # is it the weights?

    sin2theta_effective= 0.2312
    gv_ga=  -1 + 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm

    Ae=Ae_sm

    costheta_tau=math.cos(genTauP4.Theta()) # this is the theta of the Tau, not the Rho 

    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)

    Pnew_P1= - ( (+1)   * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm* (+1) *costheta_tau)
    Pnew_M1= - ( (-1)   * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm* (-1) *costheta_tau)

    ratioMass2= mtau*mtau/mRho/mRho

    term_a= 2/3*((1-Ptau_sm*z)- ratioMass2*(1+Ptau_sm*z)) + ratioMass2* (1+Ptau_sm*z)
    term_b= -2/3*((1-Ptau_sm*z-ratioMass2*(1+Ptau_sm*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Ptau_sm*math.sin(2*anglePsi)*math.sin(gen_theta_Rho))*(3*gen_cosBeta*gen_cosBeta-1)/2
    den = term_a+term_b

    term_a_P1= 2/3*((1-Pnew_P1*z)- ratioMass2*(1+Pnew_P1*z)) + ratioMass2* (1+Pnew_P1*z)
    term_b_P1= -2/3*((1-Pnew_P1*z-ratioMass2*(1+Pnew_P1*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew_P1*math.sin(2*anglePsi) *math.sin(gen_theta_Rho))*(3*gen_cosBeta*gen_cosBeta-1)/2
    num_P1=term_a_P1+term_b_P1

    term_a_M1= 2/3*((1-Pnew_M1*z)- ratioMass2*(1+Pnew_M1*z)) + ratioMass2* (1+Pnew_M1*z)
    term_b_M1= -2/3*((1-Pnew_M1*z-ratioMass2*(1+Pnew_M1*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew_M1*math.sin(2*anglePsi) *math.sin(gen_theta_Rho))*(3*gen_cosBeta*gen_cosBeta-1)/2
    num_M1=term_a_M1+term_b_M1

    weight_P1 =num_P1/den
    weight_M1 = num_M1/den

    if testAtau==0:
      return (cos_theta,cosPsi,gen_cosBeta,w,weight_P1,weight_M1)
    else:
      Pnew_test= - ( (testAtau)   * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm* testAtau  *costheta_tau)
      term_a_Ptest= 2/3*((1-Pnew_test*z)- ratioMass2*(1+Pnew_test*z)) + ratioMass2* (1+Pnew_test*z)
      term_b_Ptest= -2/3*((1-Pnew_test*z-ratioMass2*(1+Pnew_test*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew_test*math.sin(2*anglePsi) *math.sin(gen_theta_Rho))*(3*gen_cosBeta*gen_cosBeta-1)/2
      num_Ptest=term_a_Ptest+term_b_Ptest
      weight_Ptest =num_Ptest/den
      return (cos_theta,cosPsi,gen_cosBeta,w,weight_P1,weight_M1,weight_Ptest)

def wVariabRECO(RhoP4,PionP4,beamE):

    BoostRho=ROOT.TVector3()
    BoostRho=RhoP4.BoostVector()

    Pion_RhoRes=ROOT.TLorentzVector()
    Pion_RhoRes.SetPxPyPzE(PionP4.Px(),PionP4.Py(),PionP4.Pz(),PionP4.E())
    Pion_RhoRes.Boost(-BoostRho)

    mtau=1.7769
    mRho=RhoP4.M()
    mRhoRES=0.77545
    mPion=0.135

    sqrts=beamE*2
 
    x= RhoP4.E()/ beamE

#    cos_theta = 4*mtau*mtau / (mtau*mtau-mRho*mRho)* x
#    -(mtau*mtau+mRho*mRho)/(mtau*mtau-mRho*mRho) # this one does not work 

    cos_theta = (2 *x * mtau*mtau - mtau*mtau - mRho*mRho)/(mtau*mtau-mRho*mRho)/math.sqrt(1-4*mtau*mtau/sqrts/sqrts)
#    cos_theta = (2 *x * mtau*mtau - mtau*mtau - mRhoRES*mRhoRES)/(mtau*mtau-mRhoRES*mRhoRES)/math.sqrt(1-4*mtau*mtau/sqrts/sqrts)
    # careful, the good one is the invariant mass, not the pole!! do not use the commented version 

    if cos_theta>1:
     #  print ('What happened?', cos_theta)
       cos_theta=1
    if cos_theta<-1:
     #  print ('What happened?', cos_theta)
       cos_theta=-1
    angleTheta=math.acos(cos_theta)

    v1 = Pion_RhoRes.Vect()
    v2 = RhoP4.Vect()
    beta= v1.Angle(v2)
    cosBeta = math.cos( beta) # why??????

    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(abs(x*x-4*mRho*mRho/sqrts/sqrts)))
    if cosPsi>1:
     #  print ('What happened?', cosPsi)
       cosPsi=1
    if cosPsi<-1:
     #  print ('What happened?', cosPsi)
       cosPsi=-1
    anglePsi=math.acos(cosPsi)

    w_a= (-2+mtau*mtau/mRho/mRho + 2*(1+mtau*mtau/mRho/mRho)*(3*cosPsi-1)/2*(3*cosBeta*cosBeta-1)/2)* cos_theta
    w_b = 3 * math.sqrt(mtau*mtau/mRho/mRho) * (3*cosBeta*cosBeta-1)/2 * math.sin(2 * anglePsi) * math.sin(angleTheta)
    w_c = 2 + (mtau*mtau/mRho/mRho) - 2 *(1- (mtau*mtau/mRho/mRho)) * (3*cosPsi-1)/2 * (3*cosBeta*cosBeta-1)/2

    w= (w_a+w_b)/w_c

    return (cos_theta,cosPsi,cosBeta,w)


