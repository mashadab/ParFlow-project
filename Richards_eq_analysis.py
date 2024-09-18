# -*- coding: utf-8 -*-
"""
Richards equation parameters
Author:Mohammad Afzal Shadab
Email: mashadab@princeton.edu
Date: 9 Sept, 2024
"""

#import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Parameters from Celia
Ks = 0.00944 #cm/s
A  = 1.175e6
gamma = 4.74
alpha = 1.611e6#alpha = 1.611e6
theta_s = 0.287
theta_r = 0.075
gamma    = 1000#3.96

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'font.family': 'Times'})

#Colors
brown  = [181/255 , 101/255, 29/255]
red    = [255/255 ,255/255 ,255/255 ]
blue   = [ 30/255 ,144/255 , 255/255 ]
green  = [  0/255 , 166/255 ,  81/255]
orange = [247/255 , 148/255 ,  30/255]
purple = [102/255 ,  45/255 , 145/255]
brown  = [155/255 ,  118/255 ,  83/255]
tan    = [199/255 , 178/255 , 153/255]
gray   = [100/255 , 100/255 , 100/255]

#Constitutive equations
K = lambda h: Ks*A/(A + np.abs(h)**gamma)  #head in cms
C = lambda h: alpha*(theta_s - theta_r)/(alpha + np.abs(h)**gamma)**2.0*gamma*np.abs(h)**(gamma-1.0) #head in cms
D = lambda h: K(h)/C(h) #head in cms
Th= lambda h: alpha*(theta_s - theta_r)/(alpha + np.abs(h)**gamma)+theta_r #Theta, head in cms
dKdh= lambda h: Ks*gamma*A*np.abs(h)**(gamma-1.0)/(A + np.abs(h)**gamma)**2.0   #head in cms
dCdh= lambda h: gamma*alpha*(theta_s - theta_r)*(2*np.abs(h)**(gamma-1)/(alpha + np.abs(h)**gamma)**3.0 + (gamma - 1.0)*np.abs(h)**(gamma - 2.0)/(alpha + np.abs(h)**gamma)**2.0)

h = np.linspace(-100,0,100000)

'''
#plotting
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,K(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
plt.ylabel(r'$K$ [cm/s]')
plt.ylim([-0.001,0.010])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'conductivity_vs_head.pdf',bbox_inches='tight', dpi = 600)

fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,C(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
plt.ylim([-0.001,0.0065])
plt.ylabel(r'$C$ [1/cm]')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'specific_mositure_vs_head.pdf',bbox_inches='tight', dpi = 600)

fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,D(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
plt.ylim([-0.001,50])
plt.ylabel(r'$D$ [cm^2/s]')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'diffusivity_vs_head.pdf',bbox_inches='tight', dpi = 600)

fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,Th(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
#plt.ylim([-0.001,50])
plt.ylabel(r'$\theta$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'theta_vs_head.pdf',bbox_inches='tight', dpi = 600)


fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,dKdh(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
#plt.ylim([-0.001,50])
plt.ylabel(r'$dK/dh$ [1/s]')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'dKdh_vs_head.pdf',bbox_inches='tight', dpi = 600)


fig = plt.figure(figsize=(15,7.5) , dpi=100)
plt.plot(h,dCdh(h),'b-')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
#plt.ylim([-0.001,50])
plt.ylabel(r'$dC/dh$ [1/cm$^2$]')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'dCdh_vs_head.pdf',bbox_inches='tight', dpi = 600)



#Celia Sweep
#Theta
Th_haverkamp= lambda h,gamma: alpha*(theta_s - theta_r)/(alpha + np.abs(h)**gamma)+theta_r #Theta, head in cms

fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(h,Th_haverkamp(h,1),'-',color='k',alpha=0.2,label=r'$\gamma=1.0$')
plt.plot(h,Th_haverkamp(h,3.96),'-',color='k',alpha=0.4,label=r'$\gamma=3.96$')
plt.plot(h,Th_haverkamp(h,10),'-',color='k',alpha=0.6,label=r'$\gamma=10$')
plt.plot(h,Th_haverkamp(h,50),'-',color='k',alpha=0.8,label=r'$\gamma=50$')
plt.plot(h,Th_haverkamp(h,100),'-',color='k',alpha=1.0,label=r'$\gamma=100$')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$h$ [cm]')
plt.ylim([theta_r-0.001,theta_s+0.001])
plt.xlim([np.min(h),np.max(h)])
plt.ylabel(r'$\theta$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f'theta_vs_head_haverkamp.pdf',bbox_inches='tight', dpi = 600)

#Parflow - theta
alpha_vG = 1.0
Th_vG= lambda h,n: (theta_s - theta_r)/((1 + np.abs(alpha_vG*h)**n)**(1-1/n))+theta_r #Theta, head in cms

fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(Th_vG(h,1),h,'-',color='k',alpha=0.2,label=r'$n=1$')
plt.plot(Th_vG(h,2),h,'-',color='k',alpha=0.4,label=r'$n=2$')
plt.plot(Th_vG(h,10),h,'-',color='k',alpha=0.6,label=r'$n=10$')
plt.plot(Th_vG(h,50),h,'-',color='k',alpha=0.8,label=r'$n=50$')
plt.plot(Th_vG(h,100),h,'-',color='k',alpha=1.0,label=r'$n=100$')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([theta_r-0.001,theta_s+0.001])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$\theta$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f'theta_vs_head_vG.pdf',bbox_inches='tight', dpi = 600)

'''

#ParFlow computations

#Parflow - saturation
alpha_vG = 1.0
h = np.linspace(-10,0,10000)
s_s = 1.0; s_r = 0.2
sat_arr = np.linspace(s_r,s_s,10000)

#Van Genuchten
sw_vG= lambda h,n: (s_s - s_r)/((1 + np.abs(alpha_vG*h)**n)**(1-1/n))+s_r #sw, head in cms
kr_vG= lambda h,n: (1-(np.abs(alpha_vG*h)**(n-1))/(1+np.abs(alpha_vG*h)**n)**(1-1/n))**2.0/((1 + np.abs(alpha_vG*h)**n)**((1-1/n)/2)) #kr, head in cms

#Haverkamp
sw_haverkamp= lambda h,gamma: A*(s_s - s_r)/(A + np.abs(h)**gamma)+s_r #sw, head in cms
kr_haverkamp= lambda h,gamma: A/(A + np.abs(h)**gamma) #kr, head in cms

#Ideal/kinematic
sw_vG_ideal= lambda h: s_s*(h>=0)+s_r*(h<0) #sw, head in cms
kr_vG_ideal= lambda h: 1.0*(h>=0)+0.0*(h<0)
se_vG= lambda sw: ((sw - s_r)/(s_s - s_r))
kr_vG_ideal_sat = lambda sw,n: se_vG(sw)**(1/2)*(1-(1-se_vG(sw)**(1/(1-1/n)))**(1-1/n))**2
kr_haverkamp_ideal_sat = lambda sw: se_vG(sw)

fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(sw_vG(h,1),h,'-',color='k',alpha=0.2,label=r'$n=1$')
plt.plot(sw_vG(h,2),h,'-',color='k',alpha=0.4,label=r'$n=2$')
plt.plot(sw_vG(h,10),h,'-',color='k',alpha=0.6,label=r'$n=10$')
plt.plot(sw_vG(h,50),h,'-',color='k',alpha=0.8,label=r'$n=50$')
plt.plot(sw_vG(h,100),h,'-',color='k',alpha=1.0,label=r'$n=100$')
plt.plot(sw_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$s_w$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_head_vG.pdf',bbox_inches='tight', dpi = 600)

fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(kr_vG(h,1),h,'-',color='k',alpha=0.2,label=r'$n=1$')
plt.plot(kr_vG(h,2),h,'-',color='k',alpha=0.4,label=r'$n=2$')
plt.plot(kr_vG(h,10),h,'-',color='k',alpha=0.6,label=r'$n=10$')
plt.plot(kr_vG(h,50),h,'-',color='k',alpha=0.8,label=r'$n=50$')
plt.plot(kr_vG(h,100),h,'-',color='k',alpha=1.0,label=r'$n=100$')
plt.plot(kr_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([-0.01,1.01])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$k_r$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f'kr_vs_head_vG.pdf',bbox_inches='tight', dpi = 600)


fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(sw_haverkamp(h,1),h,'-',color='k',alpha=0.2,label=r'$\gamma=1$')
plt.plot(sw_haverkamp(h,2),h,'-',color='k',alpha=0.4,label=r'$\gamma=2$')
plt.plot(sw_haverkamp(h,10),h,'-',color='k',alpha=0.6,label=r'$\gamma=10$')
plt.plot(sw_haverkamp(h,50),h,'-',color='k',alpha=0.8,label=r'$\gamma=50$')
plt.plot(sw_haverkamp(h,100),h,'-',color='k',alpha=1.0,label=r'$\gamma=100$')
plt.plot(sw_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$s_w$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_head_haverkamp.pdf',bbox_inches='tight', dpi = 600)

fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(kr_haverkamp(h,1),h,'-',color='k',alpha=0.2,label=r'$\gamma=1$')
plt.plot(kr_haverkamp(h,2),h,'-',color='k',alpha=0.4,label=r'$\gamma=2$')
plt.plot(kr_haverkamp(h,10),h,'-',color='k',alpha=0.6,label=r'$\gamma=10$')
plt.plot(kr_haverkamp(h,50),h,'-',color='k',alpha=0.8,label=r'$\gamma=50$')
plt.plot(kr_haverkamp(h,100),h,'-',color='k',alpha=1.0,label=r'$\gamma=100$')
plt.plot(kr_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlabel(r'$k_r$')
plt.xlim([-0.01,1.01])
plt.ylim([np.min(h),np.max(h)])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f'kr_vs_head_haverkamp.pdf',bbox_inches='tight', dpi = 600)


kr_vG_ideal_sat = lambda sw,n: se_vG(sw)**(1/2)*(1-(1-se_vG(sw)**(1/(1-1/n)))**(1-1/n))**2
kr_haverkamp_ideal_sat = lambda sw: se_vG(sw)


fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(sw_haverkamp(h,1),kr_haverkamp(h,1),'-',color='k',alpha=0.2,label=r'$\gamma=1$')
plt.plot(sw_haverkamp(h,2),kr_haverkamp(h,2),'-',color='k',alpha=0.4,label=r'$\gamma=2$')
plt.plot(sw_haverkamp(h,10),kr_haverkamp(h,10),'-',color='k',alpha=0.6,label=r'$\gamma=10$')
plt.plot(sw_haverkamp(h,50),kr_haverkamp(h,50),'-',color='k',alpha=0.8,label=r'$\gamma=50$')
plt.plot(sw_haverkamp(h,100),kr_haverkamp(h,100),'-',color='k',alpha=1.0,label=r'$\gamma=100$')
plt.plot(sat_arr,kr_haverkamp_ideal_sat(sat_arr),'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$k_r(h)$')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([-0.01,1.01])
plt.xlabel(r'$s_w(h)$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_kr_haverkamp.pdf',bbox_inches='tight', dpi = 600)


fig = plt.figure(figsize=(7.5,5) , dpi=100)
#plt.plot(sw_vG(h,1),kr_vG(h,1),'-',color='k',alpha=0.2,label=r'$n=1$')
plt.plot(sw_vG(h,2),kr_vG(h,2),'-',color='k',alpha=0.4,label=r'$n=2$')
plt.plot(sw_vG(h,10),kr_vG(h,10),'-',color='k',alpha=0.6,label=r'$n=10$')
plt.plot(sw_vG(h,50),kr_vG(h,50),'-',color='k',alpha=0.8,label=r'$n=50$')
plt.plot(sw_vG(h,100),kr_vG(h,100),'-',color='k',alpha=1.0,label=r'$n=100$')
#plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,1),'--',color='r',alpha=0.2,label=r'kinematic')
plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,2),'--',color='r',alpha=0.4,label=r'kinematic')
plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,10),'--',color='r',alpha=0.6)
plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,50),'--',color='r',alpha=0.8)
plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,100),'--',color='r',alpha=1.0)
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$k_r(h)$')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([-0.01,1.01])
plt.xlabel(r'$s_w(h)$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_kr_vG.pdf',bbox_inches='tight', dpi = 600)


N_saturation = 100.0; N_RelPerm = 100.0
sw_vG= lambda h,n: (s_s - s_r)/((1 + np.abs(alpha_vG*h)**n)**(1-1/n))+s_r
sw_vG_invert= lambda h: sw_vG(h,N_saturation)-1.0
x = fsolve(sw_vG_invert,0.5)

kr_vG= lambda h,n: (1-(np.abs(alpha_vG*h)**(n-1))/(1+np.abs(alpha_vG*h)**n)**(1-1/n))**2.0/((1 + np.abs(alpha_vG*h)**n)**((1-1/n)/2)) #kr, head in cms
print('kr=',kr_vG(x,N_RelPerm))

print('head',x,' m')


#Derivative of kr with se
#van-Genuchten
n = 10000
m = 1 - 1/n
init_sat = 0.765685424949238
kr_vG= lambda h,n: (1-(np.abs(alpha_vG*h)**(n-1))/(1+np.abs(alpha_vG*h)**n)**(1-1/n))**2.0/((1 + np.abs(alpha_vG*h)**n)**((1-1/n)/2)) #kr, head in cms
sw_vG= lambda h,n: (s_s - s_r)/((1 + np.abs(alpha_vG*h)**n)**(1-1/n))+s_r
se_vG= lambda sw: ((sw - s_r)/(s_s - s_r))
kr_sat_vG= lambda sw: se_vG(sw)**(1/2)*(1-(1-se_vG(sw)**(1/m))**m)**2

#derivatives
delta_sw = 1e-10
dkrbydsw_vG_numeric = lambda sw: ( kr_sat_vG(sw +delta_sw) - kr_sat_vG(sw))/(delta_sw)  #numerical
dkrbydsw_vG_analytic= lambda sw: 1/(s_s - s_r) * (2*se_vG(sw)*((se_vG(sw))**(1/m)*(1 - (se_vG(sw))**(1/m))**m + ((se_vG(sw))**(1/m) - 1)*((1 - (se_vG(sw))**(1/m))**m - 1))*((1 - (se_vG(sw))**(1/m))**m - 1)/((se_vG(sw))**(1/m) - 1))  #analytical
dkrbydsw_vG_analytic2= lambda sw: 1/(s_s - s_r) * (2*se_vG(sw)**(1/m)*(1 - se_vG(sw)**(1/m))**m*(1 - (1 - se_vG(sw)**(1/m))**m)/(se_vG(sw)**0.5*(1 - se_vG(sw)**(1/m))) + 0.5*(1 - (1 - se_vG(sw)**(1/m))**m)**2/se_vG(sw)**0.5)

#Simplified
dkrbydsw_vG_analytic2simplified= lambda sw: 1/(s_s - s_r) * (0.5*(se_vG(sw))**0.5*((se_vG(sw))**(1/m) - 1)*((1 - (se_vG(sw))**(1/m))**m - 1) + 2*(se_vG(sw))**((0.5*m + 1)/m)*(1 - (se_vG(sw))**(1/m))**m)*((1 - (se_vG(sw))**(1/m))**m - 1)/((se_vG(sw))**1.0*((se_vG(sw))**(1/m) - 1))

#Non-simplified Analytic 3
dkrbydsw_vG_analytic3= lambda sw: 1/(s_s - s_r) * ((-(se_vG(sw))*(1 - (se_vG(sw))**(1/m))**m + (se_vG(sw)))*(2*(se_vG(sw))**(1/m)*(1 - (se_vG(sw))**(1/m))**m/(1 - (se_vG(sw))**(1/m)) - 2*(1 - (se_vG(sw))**(1/m))**m + 2))

#Simplified Analytic 3
dkrbydsw_vG_analytic3simplified= lambda sw: 1/(s_s - s_r) * (2*(se_vG(sw))*((se_vG(sw))**(1/m)*(1 - (se_vG(sw))**(1/m))**m + ((se_vG(sw))**(1/m) - 1)*((1 - (se_vG(sw))**(1/m))**m - 1))*((1 - (se_vG(sw))**(1/m))**m - 1)/((se_vG(sw))**(1/m) - 1))

    #1/(s_s - s_r) * (2* (se_vG(sw) - (se_vG(sw))*(1 - (se_vG(sw))**(1/m))**m) * (1 - (1-(se_vG(sw))**(1/m))**m + (se_vG(sw))**(1/m)*(1-(se_vG(sw))**(1/m))) )


fig = plt.figure(figsize=(5,7.5) , dpi=100)
#sat_array = np.linspace(s_r,s_s,10000)  #Saturation array
sat_array = np.linspace(s_r,init_sat,10000)  #Saturation array
plt.plot(sat_array, dkrbydsw_vG_numeric(sat_array),color='k',linestyle='-',label='Numeric',alpha=0.5)   
plt.plot(sat_array, dkrbydsw_vG_analytic(sat_array),'r--',label='Analytic pen-paper')
plt.plot(sat_array, dkrbydsw_vG_analytic2(sat_array),'c-',label='Analytic2')
plt.plot(sat_array, dkrbydsw_vG_analytic2simplified(sat_array),'b--',label='Analytic2_simplified')
plt.plot(sat_array, dkrbydsw_vG_analytic3(sat_array),'m-',label='Analytic3')
plt.plot(sat_array, dkrbydsw_vG_analytic3simplified(sat_array),'g-.',label='Analytic3 simplified')
plt.ylabel(r'd$k_r/$d$s_w$')
plt.xlabel(r'$s_w$')
plt.xlim([s_r,s_s])
#plt.ylim([-0.01,1.01])
plt.tight_layout()
plt.legend(loc='best',frameon=False)
plt.savefig(f's_w_vs_dkrdsw_vG.pdf',bbox_inches='tight', dpi = 600)


print(np.linalg.norm(dkrbydsw_vG_numeric(sat_array)-dkrbydsw_vG_analytic(sat_array)))


print('Analytic 2and 2 simplified work match with the numeric results (work) but Analytic 3 and 3 simplified match with pen and paper solution. ')



#Effect of alpha on the simulations
h = np.linspace(-10,0,10000)
s_s = 1.0; s_r = 0.2
sat_arr = np.linspace(s_r,s_s,10000)

#Van Genuchten
sw_vG_alpha= lambda h,n,alpha: (s_s - s_r)/((1 + np.abs(alpha*h)**n)**(1-1/n))+s_r #sw, head in cms
kr_vG_alpha= lambda h,n,alpha: (1-(np.abs(alpha*h)**(n-1))/(1+np.abs(alpha*h)**n)**(1-1/n))**2.0/((1 + np.abs(alpha*h)**n)**((1-1/n)/2)) #kr, head in cms


fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(sw_vG_alpha(h,2,1),h,'-',color='k',alpha=0.2,label=r'$\alpha=1$')
plt.plot(sw_vG_alpha(h,2,10),h,'-',color='k',alpha=0.4,label=r'$\alpha=10$')
plt.plot(sw_vG_alpha(h,2,100),h,'-',color='k',alpha=0.6,label=r'$\alpha=100$')
plt.plot(sw_vG_alpha(h,2,1000),h,'-',color='k',alpha=0.8,label=r'$\alpha=1000$')
plt.plot(sw_vG_alpha(h,2,10000),h,'-',color='k',alpha=1.0,label=r'$\alpha=10000$')
plt.plot(sw_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$s_w$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_head_vG_different_alpha_N=2.pdf',bbox_inches='tight', dpi = 600)


fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(kr_vG_alpha(h,2,1),h,'-',color='k',alpha=0.2,label=r'$\alpha=1$')
plt.plot(kr_vG_alpha(h,2,10),h,'-',color='k',alpha=0.4,label=r'$\alpha=10$')
plt.plot(kr_vG_alpha(h,2,100),h,'-',color='k',alpha=0.6,label=r'$\alpha=100$')
plt.plot(kr_vG_alpha(h,2,1000),h,'-',color='k',alpha=0.8,label=r'$\alpha=1000$')
plt.plot(kr_vG_alpha(h,2,10000),h,'-',color='k',alpha=1.0,label=r'$\alpha=10000$')
plt.plot(kr_vG_ideal(h),h,'--',color='r',alpha=1.0,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$h$ [cm]')
plt.xlim([-0.01,1.01])
plt.ylim([np.min(h),np.max(h)])
plt.xlabel(r'$k_r$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f'kr_vs_head_vG_different_alpha_N=2.pdf',bbox_inches='tight', dpi = 600)


sat_arr = np.linspace(s_r,s_s,100000)
h = np.linspace(-100,0,10000000)
fig = plt.figure(figsize=(7.5,5) , dpi=100)
plt.plot(sw_vG_alpha(h,2,1),kr_vG_alpha(h,2,1),'-',color='k',alpha=0.2,label=r'$\alpha=1$')
plt.plot(sw_vG_alpha(h,2,10),kr_vG_alpha(h,2,10),'-',color='k',alpha=0.4,label=r'$\alpha=10$')
plt.plot(sw_vG_alpha(h,2,100),kr_vG_alpha(h,2,100),'-',color='k',alpha=0.6,label=r'$\alpha=100$')
plt.plot(sw_vG_alpha(h,2,1000),kr_vG_alpha(h,2,1000),'X',color='k',alpha=0.8,label=r'$\alpha=1000$')
plt.plot(sw_vG_alpha(h,2,10000),kr_vG_alpha(h,2,10000),'o',color='k',alpha=1.0,label=r'$\alpha=10000$')
#plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,1),'--',color='r',alpha=0.2,label=r'kinematic')
plt.plot(sat_arr,kr_vG_ideal_sat(sat_arr,2),'--',color='r',alpha=1,label=r'kinematic')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$k_r(h)$')
plt.xlim([s_r-0.01,s_s+0.01])
plt.ylim([-0.01,1.01])
plt.xlabel(r'$s_w(h)$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.legend()
plt.savefig(f's_w_vs_kr_vG_different_alpha_N2.pdf',bbox_inches='tight', dpi = 600)