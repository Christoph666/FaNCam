import numpy as np

def sigma_to_x(sigma, sigma_err):
    x12 = 1. / sigma
    x12_err = x12 * sigma_err / sigma
    return (x12, x12_err)

def x12_to_sigma(x12, x12_err):
    sigma = 1. / (x12  / np.log(2))
    sigma_err = sigma * x12_err / x12
    return (sigma, sigma_err)

def distance_in_sigma(sigma1, sigma2, sigma1_err, sigma2_err):
    return abs(sigma1 - sigma2) / np.sqrt(sigma1_err**2 + sigma2_err**2)

def relative_distance(sigma1, sigma2):
    return abs(sigma1 - sigma2) / sigma1

# wysotki values
sigmaAl=0.1233
sigmaFe=0.2203
sigmaCo=0.2412
sigmaH2O=0.1097
sigmaSand=0.04012
# values computed by wysotzki
sigmaAlComp=0.155
sigmaFeComp=0.289
sigmaCoComp=0.309
sigmaH2OComp=0.235

sigmaAl_err=0.02227
sigmaFe_err=0.01567
sigmaCo_err=0.02198
sigmaH2O_err=0.02389
sigmaSand_err=0.03544
sigmaAlComp_err=np.sqrt(2.)*0.001
sigmaFeComp_err=np.sqrt(0.002**2 + 0.004**2)
sigmaCoComp_err=0.004
sigmaH2OComp_err=np.sqrt(0.001**2 + 0.006**2)

# my measured values
myx12Al=5.51
myx12Fe=3.48
myx12Co=3.49
myx12H2O=5.54
myx12Sand=7.37

myx12Al_err=0.4
myx12Fe_err=0.1
myx12Co_err=0.09
myx12H2O_err=0.72
myx12Sand_err=1.43

mysigmaAl, mysigmaAl_err = x12_to_sigma(myx12Al, myx12Al_err)
mysigmaFe, mysigmaFe_err = x12_to_sigma(myx12Fe, myx12Fe_err)
mysigmaCo, mysigmaCo_err = x12_to_sigma(myx12Co, myx12Co_err)
mysigmaH2O, mysigmaH2O_err = x12_to_sigma(myx12H2O, myx12H2O_err)
mysigmaSand, mysigmaSand_err = x12_to_sigma(myx12Sand, myx12Sand_err)

print "Al:\t\t", sigma_to_x(sigmaAl, sigmaAl_err), "\t\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaAl, sigmaAl, mysigmaAl_err, sigmaAl_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaAl, sigmaAl)), "%"
print "Fe:\t\t", sigma_to_x(sigmaFe, sigmaFe_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaFe, sigmaFe, mysigmaFe_err, sigmaFe_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaFe, sigmaFe)), "%"
print "Co:\t\t", sigma_to_x(sigmaCo, sigmaCo_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaCo, sigmaCo, mysigmaCo_err, sigmaCo_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaCo, sigmaCo)), "%"
print "H2O:\t\t", sigma_to_x(sigmaH2O, sigmaH2O_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaH2O, sigmaH2O, mysigmaH2O_err, sigmaH2O_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaH2O, sigmaH2O)), "%"
print "Sand:\t\t", sigma_to_x(sigmaSand, sigmaSand_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaSand, sigmaSand, mysigmaSand_err, sigmaSand_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaSand, sigmaSand)), "%"
print "Al comp:\t", sigma_to_x(sigmaAlComp, sigmaAlComp_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaAl, sigmaAlComp, mysigmaAl_err, sigmaAlComp_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaAl, sigmaAlComp)), "%"
print "Fe comp:\t", sigma_to_x(sigmaFeComp, sigmaFeComp_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaFe, sigmaFeComp, mysigmaFe_err, sigmaFeComp_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaFe, sigmaFeComp)), "%"
print "Co comp:\t", sigma_to_x(sigmaCoComp, sigmaCoComp_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaCo, sigmaCoComp, mysigmaCo_err, sigmaCoComp_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaCo, sigmaCoComp)), "%"
print "H2O comp:\t", sigma_to_x(sigmaH2OComp, sigmaH2OComp_err), "\tdist:\t", "{0:.2f}".format(distance_in_sigma(mysigmaH2O, sigmaH2OComp, mysigmaH2O_err, sigmaH2OComp_err)), "sig", "\trel. dist.:\t", "{0:.2f}".format(relative_distance(mysigmaH2O, sigmaH2OComp)), "%"
