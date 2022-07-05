from http.client import HTTPResponse
from telnetlib import DM
from django.http import HttpResponse
from django.shortcuts import render
import math
from math import*
from django import forms
import pandas as pd

# Create your views here.


def home(request):
    return render(request, 'home.html')

def add(request):
    '''mass of the neutral molecule entering'''
    val1 = request.GET['mol']
    enter_value = request.GET['ion']
    '''data = [['87-40-1', '2,4,6-trichloroanisole', 211.5],['607-99-8', '2,4,6-tribromoanisole', 344.8],
            ['87-86-5', 'pentachlorophenol', 266.3],['608-71-9', 'pentabromophenol', 488.6],
            ['88-06-2', '2,4,6-trichlorophenol', 197.4]]'''
    df = pd.read_csv('data.csv')
    #print(df)
    '''df = pd.DataFrame(data,columns=['CAS','element', 'mass'])'''
    df1 = df.to_string(index = False)
    df.index += 1
    dataframe = df.to_html(justify='center')
    #print(dataframe)
    '''print(df1)'''
    '''name = input("Enter CAS number: ")'''
    rslt_df = df[df['CAS'] == val1]
    rslt_df1 = rslt_df.to_string(index = False)
    #f = open("data.html", "w")

    '''print(rslt_df1)'''
    blankIndex=[''] * len(rslt_df)
    rslt_df.index=blankIndex
    #print(rslt_df)
    m = rslt_df['Mass']
    dipole_moment1 = rslt_df['Dipole Moment']
    dipole_moment = float(dipole_moment1)
    #dipole_moment = "{:.2f}".format(dipole_moment2)
    plz1 = rslt_df['Polarizability']
    plz = float(plz1)
    #plz = "{:.2f}".format(plz2)
    #print(plz)
    #f = open("/home/manjeet/project_rate/ratecal/static/data.html", "w")
    #f.write(dataframe)
    #f.close()
    
    '''mass of the ion entering'''
    mass_ion = {
        'h3o+':19.02,
        'nh4+':18.04,
        'no+': 30.01,
        'o2+':31.98
        }
    
    '''enter_value = input("Enter name of the molecule: ")'''
    for key,value in mass_ion.items():
        if key == enter_value:
            value12 = value
            '''print (value1)'''
    

    ebyn = int(request.GET['reduced'])
    temp = int(request.GET['temp'])
    '''dipole_moment = float(request.GET['dm'])
    plz = float(request.GET['plz'])'''
    alpha = plz*10**(-30)
    DM = dipole_moment*3.336*10**(-30)

    #==============================velocity from E/N====v = K_0*N_0*E/N=========
    K_o = 2.81 #cm^2 V^-1 S^-1
    N_o = 2.687*10**(19) # 10 19 cm −3
    #EoverN = float(input("\n Enter the reduced electric field (in Td): "))  #1 Td=10^(-17) V cm^2
    redu_field = ebyn*10**(-17)
    v = K_o*N_o*redu_field*(10**(-2))         #v in ms^(-1)
    #print ("drift velocity is :", v, "ms^-1")
    #================================extra conversion factors===================
    mmol = m
    mion = value12
    T = temp
    vel = v*v
    avg_n_o = 28.8
    mass_conv = 1.67*10**(-27)
    k  = 1.38*10**-23
    #=================================KE-ion=====================================
    m1 = mion*mass_conv
    M1 = avg_n_o*mass_conv
    E1 = 0.5*(m1*vel)
    E2 = 0.5*(M1*vel)
    E3 = 1.5*(k*T)
    E3_ev = E3*6.24*10**(18)            #conversion to eV
    KE_ion = E1+E2+E3
    KE_ion_ev = KE_ion*6.24*10**(18)
    KE_ion_ev1 = "{:.4f}".format(KE_ion_ev)
    #================================KE-COM======================================
    ms = (mmol)/(mmol+mion)
    ke = (KE_ion_ev-E3_ev)
    KE_cm = E3_ev+(ms*ke)
    KE_cm2 = float(KE_cm)
    KE_cm1 = "{:.4f}".format(KE_cm2)
    #==================================k-Lang====================================
    mu = (mion*mmol)/(mion+mmol)
    mu1 = mu*1.67*10**(-27)
    epsilon = 8.854*10**(-12)     # J^-1.C^2.m^-1
    q = 1.602*10**(-19)           # C
    #mu = 2.38*10**(-26)           # reduced mass
    part1 = math.pi*alpha*q**2
    part2 = mu1*epsilon
    K_L = math.sqrt(part1/part2)
    K_L1 =10**(6)*K_L

    
    res  = "{:.2f}".format(v)
    kinetic_energy_ion = "{:.2f}".format(KE_ion_ev)
    kinetic_energy_com = KE_cm
    rate_Lang = "{:.2e}".format(K_L1)
    #print(rate_Lang)
    #===============================T_effective==============================
    second_term = vel/(3*k)                                                 #average of N2 & O2
    mass = (mmol)*(mion+avg_n_o)/(mion+mmol)
    mass_1 = mass*mass_conv
    T_eff1 = T+second_term*mass_1
    T_eff = int(T_eff1)
    #==============================values-neede===========================
    c1 = 0.727143
    c2 = 3.71823
    c3 = 0.586920
    c4 = 4.97894
    tau = dipole_moment/math.sqrt(plz*T)
    eps = dipole_moment/math.sqrt(plz*KE_cm)
    #print(eps)
    theta = c3*(c4+math.log(tau))
    #==============================k_cap-final================================
    if eps > 1.5:
        S = math.exp(-2*(eps-1.5))
        a = tau**(0.4)
        b = eps*eps
        K_c1 = c1*a*b*S
        K_c2 = c2*(1-S)
        K_c3 = math.sin(theta)
        K_c4 = tau**(0.6)
        K_c5 = math.sqrt(eps-0.5)
        K_c6 = K_c2*K_c3*K_c4*K_c5
        K_c = 1+K_c1+K_c6
        K_cap = K_L1*K_c
#        print ("..............>>>>>\n      The Capture rate constant k_cap is :  ", K_cap, "10^(-9)cm^3 s^-1", "@", i, "m/s")
    elif eps <= 1.5:
        S = 1
        a = tau**(0.4)
        b = eps*eps
        K_c1 = c1*a*b*S
        K_c2 = c2*(1-S)*math.sin(theta)
        K_c = 1+K_c1+K_c2
        K_cap = K_L1*K_c
#        print ("..............>>>>>\n      The Capture rate constant k_cap is :  ", K_cap, "10^(-9)cm^3 s^-1", "@", i, "m/s")

    #print(K_cap)
    K_cap = "{:.2e}".format(K_cap)
    #print(K_cap)
    return render(request, 'results.html', {'dipole': dipole_moment, 'polar': plz, 'keion': KE_ion_ev1, 'temper': T, 'keioncom': KE_cm1, 'teff': T_eff, 'rates': K_cap, 'lang': rate_Lang})