#-- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pylab import style
# style.use('ggplot')
plt.rcParams['font.sans-serif'] = ['SimHei']
# plt.rcParams['axes.unicode_minus'] = False

'''BS EN 12354:2000强制辐射因子σf'''
def radiation_forced(f,a,b,speed_sound=340):
    k0=2*np.pi*f/speed_sound
    lambda_f=-0.964-(0.5+b/(np.pi*a))*np.log(b/a)+5*b/(2*np.pi*a)-1/(4*np.pi*a*b*k0**2)
    the_radiation_factor_forced_transmission=0.5*(np.log(k0*np.sqrt(a*b))-lambda_f)
    the_radiation_factor_forced_transmission[the_radiation_factor_forced_transmission>2]=2
    # if the_radiation_factor_forced_transmission>2 :
    #     the_radiation_factor_forced_transmission=2
    return the_radiation_factor_forced_transmission

'''BS EN 12354:2000自由弯曲波辐射因子σ'''
def radiation_free(f,a,b,thickness,density,speed_sound,Youngs_modulus, Poissons_ratio):
    fc=Critical_frequency(thickness,density,speed_sound,Youngs_modulus, Poissons_ratio)
    f11=Resonance_frequency(a,b,thickness,density,Youngs_modulus, Poissons_ratio)
    sigma=[]
    # sigma1=1/np.sqrt(1-fc/f)
    # sigma2=4*a*b*(f/speed_sound)**2
    # sigma3=np.sqrt(2*np.pi*f*(a+b)/(16*speed_sound))
    if f11<=fc/2:
        for i in f:
            if i>=fc:
                sigma1 = 1 / np.sqrt(1 - fc / i)
                sigma_p=sigma1
                # sigma.append(min(sigma1,2) )
            elif i<fc:
                lambda_1=np.sqrt(i/fc)
                delta1=((1-lambda_1**2)*np.log((1+lambda_1)/(1-lambda_1))+2*lambda_1)/(4*np.pi**2*(1-lambda_1**2)**1.5)
                if i>fc/2:
                    delta2=0
                else:
                    delta2=8*speed_sound**2*(1-2*lambda_1**2)/(speed_sound**2*np.pi**4*a*b*lambda_1*np.sqrt(1-lambda_1**2))
                sigma_p=2*(a+b)*speed_sound*delta1/(a*b*fc)+delta2
                sigma2 = 4 * a * b * (i / speed_sound) ** 2
                if (i < f11 < fc/2):
                    sigma_p=min(sigma2,sigma_p)
            sigma.append(min(sigma_p,2))
    elif f11>fc/2:
        for j in f:
            sigma1 = 1 / np.sqrt(1 - fc / f)
            sigma2=4*a*b*(f/speed_sound)**2
            sigma3=np.sqrt(2*np.pi*f*(a+b)/(16*speed_sound))
            if (j<fc) and (sigma2 <sigma3):
                sigma_p=min(sigma2,2)
            elif (j>fc) and (sigma1 <sigma3):
                sigma_p = min(sigma1, 2)
            else:
                sigma_p = min(sigma3, 2)
            sigma.append(sigma_p)
    sig=np.array(sigma)
    return sig

'''BS EN 12354:2000 实验室总耗损系数ηtot'''
def total_loss_factor(f,thickness,density,damping=0.01):
    eta_tot=damping+thickness*density/(485*np.sqrt(f))
    return eta_tot

'''拉梅弹性常数λ=lame_lambd=Eσ/[(1+σ)*(1-2σ)]，μ=lame_mju=E/[2(1+σ)]'''
def _Lame_coefficient(Youngs_modulus=2.1e11, Poissons_ratio=0.0001):
    lame_lambd = (Youngs_modulus * Poissons_ratio) / (
                (1 + Poissons_ratio) * (1 - 2 * Poissons_ratio))
    lame_mju = Youngs_modulus / ((1 + Poissons_ratio) * 2)

    return lame_lambd, lame_mju

'''纵波声速c1=sqrt((λ+2*μ)/ρp),横波声速c2=sqrt(μ/ρp)'''
def _wave_velocity(density,Youngs_modulus=2.1e11, Poissons_ratio=0.0001):

        lambd, mju=_Lame_coefficient(Youngs_modulus, Poissons_ratio)

        the_Longitudinal_wave = np.sqrt((lambd + 2 * mju) /density)
        the_Shear_wave = np.sqrt(mju / density)

        return the_Longitudinal_wave, the_Shear_wave

'''板的弯曲劲度Bs=(t^3 * ρp * Cb^2 )/12,板中的纵波速度Cb=sqrt(E/[(1-σ^2)ρp])'''
def bending_stiffness(thickness,density,Youngs_modulus=2.1e11, Poissons_ratio=0.0001 ):

        the_Plate_longitudinal_wave = np.sqrt(Youngs_modulus / ((1 -Poissons_ratio ** 2) * density))

        the_bending_stiffness=(Youngs_modulus*thickness**3)  / (12*(1-Poissons_ratio**2))

        return the_bending_stiffness,the_Plate_longitudinal_wave

"""最小临界频率frp或者fc=c^2 * sqrt(Bs/m) / (2*π)"""
def Critical_frequency(thickness,density,speed_sound,Youngs_modulus, Poissons_ratio):
    mass_area=thickness*density
    bs,cb=bending_stiffness(thickness,density,Youngs_modulus, Poissons_ratio)
    the_Critical_frequency=(speed_sound**2)*np.sqrt(mass_area/bs )/(2*np.pi)
    return the_Critical_frequency
    # return 40
"""最小共振频率fr=π*sqrt(Bs/m)*(1/a^2+1/b^2)/2=π*Cb*t*(1/a^2+1/b^2)/(4*sqrt(3))"""
def Resonance_frequency(a,b,thickness,density,Youngs_modulus, Poissons_ratio):
    bs,cb=bending_stiffness(thickness,density,Youngs_modulus, Poissons_ratio)
    fr_min=np.pi*cb* thickness*(1/a**2+1/b**2)/(4*np.sqrt(3))
    return fr_min

"""经典质量定律a=ω*m/(2*ρ*c)=π* m *f / ( ρ* c )"""
def transmission_a( f,thickness,density,speed_sound=340,air_density=1.22):
    a = f * thickness*density * np.pi / (speed_sound*air_density)
    return a

'''垂直入射质量定律R0='''
def quality_of_law(f,thickness,density,speed_sound=340,air_density=1.22):
    a = transmission_a(f,thickness,density,speed_sound,air_density)
    tau0 = 1 / ( a ** 2)
    ##角度在0度
    R0 = 10 * np.log10(1 / tau0)
    return R0

'''无规入射的质量定律'''
def random_incident(f,thickness,density,speed_sound=340,air_density=1.22):
    a= transmission_a(f,thickness,density,speed_sound,air_density)
    tau = np.log(1+a**2)/a**2
    Rr=10 * np.log10(1 / tau)
    return Rr

'''场入射的质量定律'''
def field_incident(f,thickness,density,speed_sound=340,air_density=1.22):
    Rf=quality_of_law(f,thickness,density,speed_sound,air_density)-5
    return Rf

'''苏联-扎鲍罗夫理论'''
def bzlf_insul(f,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    fc = Critical_frequency(thickness, density, speed_sound, Youngs_modulus, Poissons_ratio)
    f1 = f[f< fc]
    R1=field_incident(f1,thickness,density,speed_sound,air_density)
    f2 = f[f>=fc]
    a = 30 * np.log10(f2 / fc) - 10 * np.log10(1 / damping) - 3
    R2= 20 * np.log10(density * thickness * np.pi * fc / (air_density * speed_sound)) +a
    R = np.concatenate((R1, R2))
    return R
'''Sharp 1978理论'''
def sharp_insul(f,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    fc=Critical_frequency(thickness,density,speed_sound,Youngs_modulus, Poissons_ratio)
    if fc>=50 and fc<=8000:
        f1 = f[f < fc/2] #小于吻合频率一半的频率
        zh=f[f<fc]
        f2= zh[zh>=fc/2]
        f3=f[f > fc]#大于吻合频率的频率
        b3=10*np.log10(f3*2*damping/(np.pi*fc))
        R3=quality_of_law(f3,thickness,density,speed_sound,air_density)+b3
        # deta_f=f1/1.414 #频率带宽取Δf=1 #暂时还不清楚
        # b1=10*np.log10(3/2+np.log(2*f1/deta_f))
        R1=quality_of_law(f1,thickness,density,speed_sound,air_density)-6

        a=(R3[0]-R1[-1])/(fc-fc/2)
        b=R3[0]-a*fc
        R2 = a*f2+b

        Rsharp=np.concatenate((R1,R2,R3))
    elif fc<50:
        # f1=f2=np.array[[]]
        f3=f
        b3 = 10 * np.log10(f3 * 2 * damping / (np.pi * fc))
        R3 = quality_of_law(f3, thickness, density, speed_sound, air_density) + b3
        Rsharp=R3
    elif fc>8000:
        if fc/2 > 8000:
            f1=f
            # f3 = f2 = np.array[[]]
            b1 = 10 * np.log10(3 / 2 + np.log(2 * f1 / 0.2316))
            R1 = quality_of_law(f1, thickness, density, speed_sound, air_density) - b1
            Rsharp = R1
        elif (fc/2 < 8000) and (fc/2>4000):
            f1 = f[f < fc / 2]  # 小于吻合频率一半的频率
            f2= f[f >= fc/2]
            # f2 = zh[zh >= fc / 2]
            # f3 = f[f > fc]  # 大于吻合频率的频率
            b3 = 10 * np.log10(fc * 2 * damping / (np.pi * fc))
            R3 = quality_of_law(fc, thickness, density, speed_sound, air_density) + b3

            b1 = 10 * np.log10(3 / 2 + np.log(2 * f1 / 0.2316))
            R1 = quality_of_law(f1, thickness, density, speed_sound, air_density) - b1

            a = (R3 - R1[-1]) / (fc - fc / 2)
            b = R3 - a * fc
            R2 = a * f2 + b

            Rsharp = np.concatenate((R1, R2))
    return Rsharp

'''Brekke 1981理论'''
def brekke_insul(f,a,b,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    fc = Critical_frequency(thickness, density, speed_sound, Youngs_modulus, Poissons_ratio)
    Rad=1 #板的辐射比未知
    if fc>=50 and fc<=8000:
        f1 = f[f < fc] #小于吻合频率一半的频率
        Rr = 20 * np.log10(thickness * density) + 30 * np.log10(f1)-20 * np.log10(Rad) + 10 * np.log10(damping) - 10 * np.log10(fc) - 44
        kt=(2*np.pi*f1/speed_sound)*np.sqrt(a*b)
        Rnr = quality_of_law(f1,thickness,density,speed_sound,air_density)- 10*np.log10(np.log(kt))+20*np.log10(1-(f1/fc)**2)
        R1=Rnr
        f2 = f[f > fc]  # 大于吻合频率的频率
        R2 = 20 * np.log10(thickness * density) + 30 * np.log10(f2) -20 * np.log10(Rad) + 10 * np.log10(damping) - 10 * np.log10(fc) - 44
        Rbrekke=np.concatenate((R1,R2))
    elif fc>8000:
        kt = (2 * np.pi * f / speed_sound) * np.sqrt(a * b)
        Rnr = quality_of_law(f, thickness, density, speed_sound, air_density) - 10 * np.log10(
            np.log(kt)) + 20 * np.log10(1 - (f/ fc) ** 2)
        Rbrekke=Rnr
    else:
        Rbrekke=20 * np.log10(thickness * density) + 30 * np.log10(f) - 20 * np.log10(Rad) + 10 * np.log10(
            damping) - 10 * np.log10(fc) - 44
    return Rbrekke

'''BS EN 12354 a,σf,σ,ηtot'''
def BS_EN_factor(f,a,b,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    lambda_f=radiation_forced(f,a,b,speed_sound)
    lambda_0=radiation_free(f,a,b,thickness,density,speed_sound,Youngs_modulus, Poissons_ratio)
    # fc=Critical_frequency(thickness,density,speed_sound,Youngs_modulus, Poissons_ratio)
    eta_tot=total_loss_factor(f,thickness,density,damping)
    a= transmission_a(f, thickness, density, speed_sound, air_density)
    return lambda_f,lambda_0,eta_tot,a

'''BS EN 12354:2000 理论'''
def Bs_EN_insul(f,a,b,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):

    fc=Critical_frequency(thickness,density,speed_sound,Youngs_modulus, Poissons_ratio)
    if fc>=50 and fc<=8000:
        f1 = f[f > fc]
        lambda_f, lambda_0,eta_tot, a1=BS_EN_factor(f1,a,b,thickness,density,speed_sound,air_density,damping,Youngs_modulus, Poissons_ratio)
        tau1=np.pi*fc*lambda_0**2/(2*eta_tot*f1*a1**2)
        R1=10 * np.log10(1 / tau1)

        f2=f[f==fc]
        if f2 is None:
            R2=np.array([])
        else:
            lambda_f, lambda_0,eta_tot, a2=BS_EN_factor(f2,a,b,thickness,density,speed_sound,air_density,damping,Youngs_modulus, Poissons_ratio)
            tau2 = np.pi * lambda_0 ** 2 / (2 * eta_tot  * a2** 2)
            R2=10 * np.log10(1 / tau2)

        f3=f[f<fc]
        lambda_f, lambda_0,eta_tot, a3=BS_EN_factor(f3,a,b,thickness,density,speed_sound,air_density,damping,Youngs_modulus, Poissons_ratio)
        tau3=(2*lambda_f+(a+b)**2*np.sqrt(fc/f3)*lambda_0**2/(eta_tot*(a**2+b**2)))/a3**2
        R3=10 * np.log10(1 / tau3)
        Rbsen = np.concatenate((R3, R2, R1))
    elif fc<50:
        lambda_f, lambda_0, eta_tot, a1 = BS_EN_factor(f, a, b, thickness, density, speed_sound, air_density, damping,
                                                       Youngs_modulus, Poissons_ratio)
        tau1 = np.pi * fc * lambda_0 ** 2 / (2 * eta_tot * f * a1 ** 2)
        Rbsen = 10 * np.log10(1 / tau1)
    elif fc>8000:
        lambda_f, lambda_0, eta_tot, a1 = BS_EN_factor(f, a, b, thickness, density, speed_sound, air_density, damping,
                                                       Youngs_modulus, Poissons_ratio)
        tau3 = (2 * lambda_f + (a + b) ** 2 * np.sqrt(fc / f) * lambda_0 ** 2 / (eta_tot * (a ** 2 + b ** 2))) / a1** 2
        Rbsen = 10 * np.log10(1 / tau3)

    return Rbsen

'''Moser 理论'''
def Moser_insul(f,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    fc = Critical_frequency(thickness, density, speed_sound, Youngs_modulus, Poissons_ratio)
    if fc>=50 and fc<=8000:
        # f1 = f[f > fc]
        # a1= transmission_a(f1, thickness, density, speed_sound, air_density)
        # eta_tot = total_loss_factor(f1, thickness, density, damping)
        # tau1=(1/a1**2)/(1/a1**2+eta_tot*abs(np.cos(fc))) ########?????????????
        # R1=10 * np.log10(1 / tau1)

        f2=f[f>=fc]
        if f2 is None:
            R2=np.array([])
        else:
            a1 = transmission_a(f2, thickness, density, speed_sound, air_density)
            eta_tot = total_loss_factor(f2, thickness, density, damping)
            tau2 = np.sqrt(fc/f2) / (2 * eta_tot  * a1** 2)
            R2=10 * np.log10(1 / tau2)

        f3=f[f<fc]
        a1 = transmission_a(f3, thickness, density, speed_sound, air_density)
        tau3=1/a1**2
        R3=10 * np.log10(1 / tau3)-3
        Rmoser = np.concatenate((R3, R2))
    elif fc<50:
        a1 = transmission_a(f, thickness, density, speed_sound, air_density)
        eta_tot = total_loss_factor(f, thickness, density, damping)
        tau1 = (1 / a1 ** 2) / (1 / a1 ** 2 + eta_tot * np.cos(fc))
        Rmoser = 10 * np.log10(1 / tau1)
    elif fc>8000:
        a1 = transmission_a(f, thickness, density, speed_sound, air_density)
        tau3 = 1 / a1 ** 2
        Rmoser = 10 * np.log10(1 / tau3)-3

    return Rmoser

'''Kristensen & Rindel 理论 1989'''
def kris_factor(f,a,b, thickness, density, speed_sound, air_density,damping):
    # Rf = random_incident(f, thickness, density, speed_sound, air_density)
    # # R0 = Rf - 10 * np.log10(np.log(2 * np.pi * f * np.sqrt(a * b) / speed_sound))
    R0 = quality_of_law(f,thickness,density,speed_sound,air_density)
    eta_tot = total_loss_factor(f, thickness, density, damping)

    return R0, eta_tot

def Kris_insul(f,a,b,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
    fc = Critical_frequency(thickness, density, speed_sound, Youngs_modulus, Poissons_ratio)
    fh=(speed_sound/(6*thickness))**2/fc
    f11=Resonance_frequency(a, b, thickness, density, Youngs_modulus, Poissons_ratio)

    f1 = f[f > fh]
    R0, eta_tot=kris_factor(f1,a,b, thickness, density, speed_sound, air_density,damping)
    deta_rh=10*np.log10(np.sqrt(f1/(5*fh)+1)+(f1/(5*fh)))
    R1=R0+10*np.log10(eta_tot*f1/fc)-deta_rh-2

    f21 = f[f >= fc];    f2 = f21[f21 < fh]
    R0, eta_tot = kris_factor(f2, a, b, thickness, density, speed_sound, air_density,damping)
    R2= R0 + 10 * np.log10(eta_tot) + 10 * np.log10(f2/fc) - 2

    f31 = f[f < fc];    f3 = f31[f31 > f11]
    R0, eta_tot = kris_factor(f3, a,b, thickness, density, speed_sound, air_density,damping)
    deta_d=(0.2+np.log(2*np.pi*f3*np.sqrt(a*b)/speed_sound))/2
    R3 = R0 -10 * np.log10(2*deta_d)+20*np.log10(1-(f3/fc)**2)

    f4 = f[f < f11]
    R0, eta_tot = kris_factor(f4, a,b, thickness, density, speed_sound, air_density,damping)
    deta_d = (0.2 + np.log(2 * np.pi * f4 * np.sqrt(a * b) / speed_sound)) / 2
    R4 = R0 - 10 * np.log10(2 * deta_d)+40*np.log10(f11/f4)

    Rkris = np.concatenate((R4,R3, R2, R1))
    return Rkris

"""Davy 模型 2009"""
class Davy_insul():
    def __init__(self,f,a,b,thickness,density,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=2.1e11, Poissons_ratio=0.25):
        self.n=2
        self.w=1.3
        self.bt=0.124
        self.s=a*b
        self.e=a*b/(a+b)
        self.fc=Critical_frequency(thickness, density, speed_sound, Youngs_modulus, Poissons_ratio)
        # self.a1=transmission_a(f1, thickness, density, speed_sound, air_density)
        # self.r=f/fc

        self.R=self.RDavy(f, thickness, density, speed_sound, air_density, Poissons_ratio,  Youngs_modulus, damping)

    def k_number(self,f,speed_sound):
        k=2*np.pi*f/speed_sound
        # self.k=k
        return k

    def p_h_af_q_number(self,f,speed_sound):
        k=self.k_number(f,speed_sound)
        p=self.w*np.sqrt(np.pi/(2*k*self.e))
        p[p>1]=1
        # p=min(p1,1)
        self.p = p

        h = 1 / (2 * np.sqrt(2 * k * self.e / np.pi) / 3 - self.bt)
        self.h=h

        af=h/p-1
        # self.af=af

        q=2*np.pi/((k**2) * self.s)
        self.q=q
        return p,h,af,q

    def k_l_b_s_X_GN_G_kt_number(self,f,thickness,density,Poissons_ratio,Youngs_modulus):
        G = Youngs_modulus / (2 * (1 + Poissons_ratio))
        kl_2=((2*np.pi*f)**2) * density*(1-Poissons_ratio**2)/Youngs_modulus
        kb_4=12*kl_2/(thickness**2)
        X=((1+Poissons_ratio)/(0.87+1.12*Poissons_ratio))**2
        GN=G/X
        kt_2=((2*np.pi*f)**2)*density/GN
        ks_2=kt_2+kl_2

        # C=np.sqrt(1-kt_2*thickness**2/12+ks_2**2/(4*kb_4))
        # self.C=C
        return kl_2,kt_2,ks_2,kb_4,X

    def _g_number(self,f):
        if (np.max(f)>self.fc) and (np.min(f)<self.fc):
            raise TypeError('g频率出错')
        elif np.min(f)>= self.fc:
            g=np.sqrt(1-self.fc/f)
        elif np.max(f)<self.fc:
            nu=len(list(f))
            g=np.zeros(nu)

        return g

    def delt_z(self,f,speed_sound):
        g=self._g_number(f)
        p, h, af, q= self.p_h_af_q_number(f,speed_sound)
        z=[]
        for i in range(len(list(f))):
            if g[i]>=p[i] and g[i]<=1:
                # zt=1/np.power(g[i]**self.n+q[i]*self.n,self.n)
                zt = 1 / np.sqrt(g[i] ** 2 + q[i] ** 2)
                z.append(zt)
            elif g[i]>=0 and g[i]<p[i]:
                # zt=1/np.power((h[i]-af[i]*g[i])**self.n+q[i]*self.n,self.n)
                zt=1/np.sqrt((h[i]-af[i]*g[i])**2+q[i]**2)
                z.append(zt)
        del_z=np.array(z)

        return del_z

    def km_number(self,f,thickness,density,Poissons_ratio,Youngs_modulus):
        a=Youngs_modulus*thickness**2/((1-Poissons_ratio**2)*density*12)
        kl_2, kt_2, ks_2, kb_4, X=self.k_l_b_s_X_GN_G_kt_number(f,thickness,density,Poissons_ratio,Youngs_modulus)
        c = -(2 * np.pi * f) ** 2
        b=c*(1+2*X/(1-Poissons_ratio))*thickness**2/12
        kmt=[]
        for i in range(len(list(f))):
            if a == 0:
                # raise TypeError('a不能为0')
                km=-c[i]/b[i]
            else:
                if not isinstance(a, (int, float)) or not isinstance(b[i], (int, float)) or not isinstance(c[i], (int, float)):
                    # km=-1
                    raise TypeError('Bad operand type')
                delta1 = b[i] ** 2 - 4 * a * c[i]
                if delta1 < 0:
                    # return '无实根'
                    km=-1
                else:
                    x1 = (np.sqrt(delta1) - b[i]) / (2 * a)
                    # x1=max(x1,0)
                    x2 = -(np.sqrt(delta1) + b[i]) / (2 * a)
                    # x2=max(x2,0)
                    if x1==x2 and x1>0:
                        km=x1
                    elif x1==x2 and x1==0:
                        km=0
                    elif max(x1,x2)<0:
                        km=-1
                    # elif min(x1,x2)>0:
                    #     km=x1
                    else:
                        km=max(x1,x2)
            if km<0:
                raise TypeError('Bad operand type')
            else:
                kmt.append(km)
        return kmt

    def A_B_C_x_number(self,f,thickness,density, speed_sound, air_density,Poissons_ratio,Youngs_modulus,damping):
        kl_2, kt_2, ks_2, kb_4, X=self.k_l_b_s_X_GN_G_kt_number(f,thickness,density,Poissons_ratio,Youngs_modulus)
        del_z=self.delt_z(f,speed_sound)
        km=self.km_number(f,thickness,density,Poissons_ratio,Youngs_modulus)
        A=(1+thickness**2*(km*kt_2/kl_2-kt_2)/12)**2
        B=1-kt_2*thickness**2/12+km*ks_2/kb_4
        C = np.sqrt(1 - kt_2 * thickness ** 2 / 12 + ks_2 ** 2 / (4 * kb_4))
        a1=transmission_a(f, thickness, density, speed_sound, air_density)
        damping_t = total_loss_factor(f,thickness,density,damping)
        x=del_z+a1*damping_t
        r=f/self.fc
        return A,B,C,x,r,a1

    def RDavy(self,f, thickness, density, speed_sound, air_density, Poissons_ratio,
                                                Youngs_modulus, damping):

        f1 = f[f >=self.fc]
        A, B, C, x, r, a1 = self.A_B_C_x_number(f1, thickness, density, speed_sound, air_density, Poissons_ratio,
                                                Youngs_modulus, damping)
        del_z = self.delt_z(f1, speed_sound)
        # td=A*(del_z)**2*(np.arctan(2*a1/x)-np.arctan(2*a1*(1-r)/x))/(B*C*2*a1*r*x)
        td =  (del_z) ** 2 * (np.arctan(2 * a1 / x) - np.arctan(2 * a1 * (1 - r) / x)) / ( 2 * a1 * r * x)
        # damping_t = total_loss_factor(f1, thickness, density, damping)
        # td=np.pi/(a1**2*2*(r-1)*damping_t)
        xz=A/(B*C)
        td=td*xz
        R1=10*np.log10(1/td)

        f2 = f[f < self.fc]
        p, h, af, q=self.p_h_af_q_number(f2,speed_sound)
        A, B, C, x, r, a1 = self.A_B_C_x_number(f2, thickness, density, speed_sound, air_density, Poissons_ratio,
                                                Youngs_modulus, damping)
        del_z = self.delt_z(f2, speed_sound)
        td = (del_z) ** 2 * (np.arctan(2 * a1 / x) - np.arctan(2 * a1 * (1 - r) / x)) / (2 * a1 * r * x)
        xz = A / (B * C)
        td = td * xz
        npc=np.log((1+np.sqrt(1+q**2))/(p+np.sqrt(p**2+q**2)))+np.log((h+np.sqrt(h**2+q**2)/(p+np.sqrt(p**2+q**2))))/af
        td2=2*npc/a1**2
        td2=td+td2
        R2 = 10 * np.log10(1 / td2)
        Rdavy = np.concatenate((R2, R1))
        return Rdavy






ce=np.array([22.0, 22.0, 23.0,
             24.24,	25.99,	31.07,	29.72,	31.92	,33.69,	33.73	,35.05	,37.48	,39.73,	40.78	,42.96,	44.49,	44.10,	37.23,	37.25	,42.80,	46.11,
             47,48])
frequency = np.array([50, 63, 80,
                                   100,  125,  160,  200,  250,  315,  400,  500,  630,  800,
                                   1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000])

# q=radiation_free(frequency,1.9,2.2,0.005,7800,340,2.1e11, 0.25)

a=quality_of_law(frequency,0.005,7800)
b=random_incident(frequency,0.005,7800)
c=field_incident(frequency,0.005,7800)
cc=bzlf_insul(frequency,0.005,7800)
d=sharp_insul(frequency,0.005,7800)
e=brekke_insul(frequency,1.9,2.2,0.005,7800)
g=Bs_EN_insul(frequency,1.9,2.2,0.005,7800,damping=0.03)
# g=Bs_EN_insul(frequency,2.8,3.6,0.015,787,damping=0.09,Youngs_modulus=2.98e9, Poissons_ratio=0.22)
# g=Bs_EN_insul(frequency,2.8,3.6,0.06,1011,damping=0.03,Youngs_modulus=11.85e9, Poissons_ratio=0.15)
h=Moser_insul(frequency,0.005,7800)
l=Kris_insul(frequency,1.9,2.2,0.005,7800)
O=Davy_insul(frequency,1.9,2.2,0.005,7800,damping=0.03).R
# # O=Davy_insul(frequency,4.08,3.08,0.2,2300,speed_sound=340,air_density=1.22,damping=0.01,Youngs_modulus=1.6e10, Poissons_ratio=0.3).R

# O=Davy_insul(frequency,2.8,3.6,0.015,787,damping=0.09,Youngs_modulus=2.98e9, Poissons_ratio=0.22).R
#gypsum plaster board 13mm
# O=Davy_insul(frequency,2.44,3.05,0.013,770,speed_sound=340,air_density=1.22,damping=0.04,Youngs_modulus=0.185e10, Poissons_ratio=0.3).R


f1 =['50', '63', '80','100','125','160','200','250','315','400','500','630','800','1K','1.25K','1.6K','2K','2.5K','3.15K','4K','5K','6.3k','8k']
f=np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])

plt.figure(figsize=(9,6))#固定纵横坐标比例
plt.plot(f,ce, label='钢5mm')
# plt.plot(f,a, label="垂直")
# plt.plot(f,b, label="无规")
# plt.plot(f,c, label="场射")
plt.plot(f,cc, label="bzlf")
plt.plot(f,d, label="Sharp")
plt.plot(f,e, label="Brekke")
plt.plot(f,g, label="BS EN 12354")
plt.plot(f,h, label="Moser")
plt.plot(f,l, label="Kris")
plt.plot(f,O, label="Davy")

plt.xticks(f, f1,rotation=60)#设置横坐标，及倾斜度数
from matplotlib.font_manager import FontProperties
font = FontProperties(fname=r"C:\Windows\Fonts\simhei.ttf", size=14)
plt.title(u'隔声模拟曲线', fontproperties=font)
plt.legend()
plt.show()