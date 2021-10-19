#For QscR-Promotor and GFP
QscR_param = {
'basexpr' : 0.2, #TODO fix
'ktransl' : 2.04, #https://www.denovodna.com:4433/shared/3qEKCv9GJEMkSPjyRiyWOSrU6pw2nsIM
'ktranskr' : 2.04, # https://www.denovodna.com:4433/shared/dvXq9kDxF9aDKfUh7hxuPZoZldJ89EQz
'Cn' : 17,
'dmrna' : 0.3301,
'dprot' : 0.015, #TODO fix
'Kd_TF' : 1.213433, # https://journals.asm.org/doi/10.1128/JB.01041-10#:~:text=There%20are%20two%20acyl%2Dhomoserine,sensing%20signals%20produced%20by%20P.&text=The%20receptors%20for%20these%20quorum,9%2C%2010%2C%2035)
'Kd_TFP' : 0.861146, #0.53,
'nh' : 1,
}

#For LasR and GFP
LasR_param = {
'basexpr' : 1, #TODO fix
'ktransl' : 1.35,# https://www.denovodna.com:4433/shared/hrrBbTW49bId1dfDi0BiTSQrWeaKK3Mc
'ktranskr' : 506.38,#https://www.denovodna.com:4433/shared/f5R1axHGySFJpzlyxx1TaQGqczL8hg7Z
'Cn' : 1, #nochmal drüberschauen für E.Coli
'dmrna' : 7.197,
'dprot' : 1, #TODO fix
'Kd_TF' : 1, #TODO fix
'Kd_TFP' : 1, #TODO fix
'nh' : 1, #TODO fix
}

# For diffusion through cell membrane of E.coli
Membrane_Diff_param = {
'D': 2, # in min^-1, https://www.biorxiv.org/content/10.1101/106229v1
'Vcell' : 5.31*10**(-12),#in L, http://book.bionumbers.org/what-is-the-concentration-of-bacterial-cells-in-a-saturated-culture/
'Venv' : 1398.72*10**(-12),#in L, http://book.bionumbers.org/what-is-the-concentration-of-bacterial-cells-in-a-saturated-culture/
}

#for AHL-Regulator with Langmuir_hill
Ahl_reg_param = {
'Kd' : 1, #TODO fix
'nh' : 1, #TODO fix
}
