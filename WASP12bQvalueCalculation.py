# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 00:35:15 2020

@author: ar-h1
"""

import numpy as np
from astropy import units as u
from astropy.constants import G

Pdot = -32.53*(u.ms/u.year)

#Mp = 1.470*u.Mjup
Mp = 1.46*u.Mjup
Ms = 1.434*u.Msun
Rs = 1.657*u.Rsun
a = 0.02340*u.au

P = 1.091420090*u.day


Factor1 = -27*np.pi/(2*Pdot)

Factor2 = Mp/Ms

RsOvera = 1/3.1089

#Factor3 = (Rs/a)**5
Factor3 = (RsOvera)**5

Qprime_s = (Factor1*Factor2*Factor3).decompose()

print(r'$Q\'_\star$ = %.3e'%(Qprime_s))


Ef1 = (((2*np.pi)**(2/3))/3)*Mp

Ef2 = (G*Ms/P)**(2/3)

Ef3 = (1/P)*(Pdot)

dEdt = (Ef1*Ef2*Ef3).to(u.watt)

Lf1 = Mp/(3*(2*np.pi)**(1/3))
Lf2 = (G*Ms/P)**(2/3)
Lf3 = Pdot 

DLdt = (Lf1*Lf2*Lf3).decompose()

Vorb = ((2*np.pi*a)/P).decompose()

Eorb = (0.5*Mp*(Vorb**2)).to(u.J)

Lorb = (Mp*Vorb*a).decompose()

Etime = (Eorb/dEdt).to(u.Myr)
Ltime = (Lorb/DLdt).to(u.Myr)

DecayTimescale = (P/Pdot).to(u.Myr)


Ptide = P/2

BarkerQ = (1.5e5)*((Ms/(1*u.Msun))**2)*(Rs/(1*u.Rsun))*((Ptide/(0.5*u.day))**(8/3))
print('BarkerQ: %.3e'%(BarkerQ))

