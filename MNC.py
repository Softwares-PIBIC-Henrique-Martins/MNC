# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 10:50:39 2021

@author: Henrique Martins Ferreira
"""

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)
from matplotlib.figure import Figure
import tkinter as tk
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad_vec
import pandas as pd
from tkinter import ttk

#Funções
class A():
    
    def Variaveis(self):
        self.mt=float(eval(Entrada_mt.get()))
        self.rho=float(eval(Entrada_rho.get()))
        self.uim=float(eval(Entrada_uim.get()))
        self.g=float(eval(Entrada_g.get()))
        self.z0=float(eval(Entrada_z0.get()))
        self.z1=float(eval(Entrada_z1.get()))
        self.R=float(eval(Entrada_R.get()))
        self.Msv=float(eval(Entrada_Msv.get()))
        self.Dm=float(eval(Entrada_Dm.get()))
        self.sigmad=float(eval(Entrada_sigmad.get()))
        self.u0=float(eval(Entrada_u0.get()))
        self.H=float(eval(Entrada_H.get()))
        self.Kb=float(eval(Entrada_Kb.get()))
        self.T=float(eval(Entrada_T.get()))
        self.t0=float(eval(Entrada_t0.get()))
        self.v0=self.mt/self.rho #m**3 volume da amostra
        self.dis=(1/self.z0**4)-(1/self.z1**4)-(self.z0**2+self.R**2)/(3*(self.z0**2+self.R**2)**3)+(self.z1**2+self.R**2)/(3*(self.z0**2+self.R**2)**3)
        self.Hk_=2/(self.u0*self.Msv) #campo de anisotropia dividido por Kef
        self.h_=self.H/self.Hk_
        self.v1=6*self.Kb*self.T
        Label_Kef_c.grid_forget()
        Label_Kef_mc.grid_forget()
        self.Progresso_mc=ttk.Progressbar(frame_saida, variable=var_barra, maximum=float(Entrada_num_ensaios_mc.get()),length=300)
        self.Progresso_mc.grid_forget()
        self.Bifurc()
        
    def MonteCarlo(self,f,xs,ys,p0,raio,num_ensaios):
        Label_Kef_mc.grid(row=1, column=1)
        Label_Kef_mc["text"]="\n Começando a calcular o ajuste por Monte Carlo\n \n "
        janela.update()
        ListaIndexEnsaios=[]
        sigma=np.std(ys)
        n_dados=len(xs)
        
        Termo1=[(f(xs[j],*p0)-ys[j])**2/sigma**2 for j in range(n_dados)]
        Qui0=np.sum(Termo1)
        ListaQuis=[Qui0]
        Transposta=[2*raio[i]*np.random.random_sample(num_ensaios)-raio[i] for i in range(len(p0))]
        Ensaios=list(map(lambda *i: [j for j in i], *Transposta))
        index_ensaios=0
        Label_Kef_mc.grid(row=1,column=1)
        for delta in Ensaios:
            Termo=[(f(xs[j],*[p0[i]+delta[i] for i in range(len(p0))])-ys[j])**2/sigma**2 for j in range(n_dados)] #Termo é cada termo da soma que compõe o quiquadrado
            Qui1=np.sum(Termo)
            if Qui1 < Qui0:
                Qui0=Qui1
                ListaQuis.append(Qui0)
                ListaIndexEnsaios.append(index_ensaios)
            index_ensaios+=1
            self.Progresso_mc.grid(row=1,column=2,stick=tk.NW)
            var_barra.set(index_ensaios)
            Label_Kef_mc["text"]="Monte Carlo\nJá foram feitos "+str(index_ensaios)+" ensaios\n \n "
            janela.update()
        Label_Kef_mc.grid_forget()
        self.Progresso_mc.grid_forget()
        janela.update()
        if len(ListaIndexEnsaios)>0:
            i=ListaIndexEnsaios[-1]
            estimativa=[p0[k]+Ensaios[i][k] for k in range(len(p0))] #valor mínimo das médias (na lista med)=valor obtido para o parâmetro a
            valores_de_delta=[Ensaios[k] for k in ListaIndexEnsaios]
        else:
            estimativa=p0 #valor mínimo das médias (na lista med)=valor obtido para o parâmetro a
            valores_de_delta=[Ensaios[k] for k in ListaIndexEnsaios]
        return estimativa, valores_de_delta, ListaQuis
    
    def AG(self, f, xs, ys, p0, raio, limite):
        Label_Kef_ag.grid(column=2, row=1, stick=tk.W, padx=5)
        Label_Kef_ag['text']="\n Começando a calcular o ajuste por Algoritmo Genético\n \n "
        janela.update()
        sigma=np.std(ys)
        n_dados=len(xs)
        
        Termo1=[(f(xs[j],*p0)-ys[j])**2/sigma**2 for j in range(n_dados)]
        Qui0=np.sum(Termo1)
        ListaQuis=[Qui0]
        
        Ensaios_efetivos=[]
        Quantidade=0
        n_ensaio=0
        n_interação=0
        
        while Quantidade<limite:
            n_ensaio=n_ensaio+1
            Ensaio=[2*raio[i]*np.random.random_sample()-raio[i] for i in range(len(p0))]
            Termo=[(f(xs[j],*[p0[i]+Ensaio[i] for i in range(len(p0))])-ys[j])**2/sigma**2 for j in range(n_dados)] #Termo é cada termo da soma que compõe o quiquadrado
            Qui1=np.sum(Termo)
            Quantidade=n_ensaio-n_interação
            Label_Kef_ag['text']="Algoritmo Genético\nLimite: {0}/{1}\n \n ".format(Quantidade, limite)
            if Qui1<Qui0:
                janela.update()
                Qui0=Qui1
                ListaQuis.append(Qui0)
                Ensaios_efetivos.append(Ensaio)
                n_interação=n_ensaio
        
        if len(Ensaios_efetivos)>=1:        
            estimativa=[p0[i]+Ensaios_efetivos[-1][i] for i in range(len(p0))]
        else:
            estimativa=p0
        
        return estimativa, Ensaios_efetivos, ListaQuis

        
    def f(self, t, Kef):
        Integrando1=lambda D: np.exp(-(np.log(D/self.Dm))**2/(2*self.sigmad**2))
        Integrando2=lambda D: 1/np.tanh(self.u0*np.pi*self.Msv*self.H*D**3/self.v1)-self.v1/(self.u0*np.pi*self.Msv*self.H*D**3)
        Integrando3=lambda D: 1-np.exp(-t*((1-(self.h_/Kef)**2)*((1-self.h_/Kef)*np.exp(-Kef*np.pi*D**3/self.v1*(1-self.h_/Kef)**2)+(1+self.h_/Kef)*np.exp(-Kef*np.pi*D**3/self.v1*(1+self.h_/Kef)**2)))/(self.t0))
        Integrando4=lambda D: 1-np.exp(-self.xs[len(self.xs)-1]*((1-(self.h_/Kef)**2)*((1-self.h_/Kef)*np.exp(-Kef*np.pi*D**3/self.v1*(1-self.h_/Kef)**2)+(1+self.h_/Kef)*np.exp(-Kef*np.pi*D**3/self.v1*(1+self.h_/Kef)**2)))/(self.t0))
        ma=quad_vec(lambda D: Integrando1(D)*Integrando2(D)*Integrando3(D),self.Dm/10,10*self.Dm)[0]
        ma1=quad_vec(lambda D: Integrando1(D)*Integrando2(D)*Integrando4(D),self.Dm/10,10*self.Dm)[0]
        return ma/ma1
    
    def Bifurc(self):
        if var_MonteCarlo.get()==1:
            res=tk.messagebox.askyesno(title="Ajuste por Monte Carlo", message="Tem certeza que quer fazer um ajuste por Monte Carlo?\nIsso pode ser demorado.")
            if res==True:
                if var_AG.get()==1:
                    res1=tk.messagebox.askyesno(title="Ajuste por Algoritmo Gnético", message="Tem certeza que quer fazer um ajuste por Algoritmo Genético?\nIsso pode ser demorado.")
                    if res1==True:
                        self.Plot()
                else:
                    self.Plot()
        elif var_AG.get()==1:
            res1=tk.messagebox.askyesno(title="Ajuste por Algoritmo Gnético", message="Tem certeza que quer fazer um ajuste por Algoritmo Genético?\nIsso pode ser demorado.")
            if res1==True:
                self.Plot()
        else:
            self.Plot()
    
    def Plot(self):
        grafico.clear()
        linhas=dadosdig.get(1.0, "end-1c").split("\n")
        self.xs=[]
        self.ys=[]
        for j in linhas:
            if j!="":
                j_dot=j.replace(",",".")
                x, y=j_dot.split("\t")
                self.xs.append(float(x))
                self.ys.append(float(y))
        self.t=np.arange(0,self.xs[-1],10)
        self.ys=[i/(self.ys[len(self.ys)-1]) for i in self.ys]
        grafico.plot(self.xs,self.ys, 'o')
        grafico.set_xlabel("Tempo (s)")
        grafico.set_ylabel("Massa normalizada")
        canvas.draw_idle()
        self.Ajuste()
        
    def curve_fit(self):
        self.est_c, self.erro_c = curve_fit(self.f, self.xs, self.ys, p0=float(eval(Entrada_p0.get())))
        grafico.plot(self.t, self.f(self.t,self.est_c[0]), '-', label="Curve_fit")
        grafico.legend()
        canvas.draw_idle()
        Label_Kef_c.grid(column=0, row=1, stick=tk.W, padx=5)
        Label_Kef_c["text"]="Mínimos Quadrados\nKef ={0} J/m³\nQuiquadrado = ...\nErro ={1} J/m³".format(self.est_c[0], self.erro_c[0][0])
        janela.update()
        self.sigma=np.std(self.ys)
        self.Termo=[(self.f(self.xs[j],*self.est_c)-self.ys[j])**2/self.sigma**2 for j in range(len(self.xs))]
        self.Qui=np.sum(self.Termo)
        print("Kef_c =", self.est_c[0])
        Label_Kef_c["text"]="Mínimos Quadrados\nKef ={0} J/m³\nQuiquadrado ={1}\nErro ={2} J/m³".format(self.est_c[0],self.Qui, self.erro_c[0][0])
        
    def MonteCarlo_fit(self):
        self.est_mc, self.lista, self.listaquis = self.MonteCarlo(self.f, self.xs, self.ys, [float(eval(Entrada_p0_mc.get()))], [float(eval(valor_raio_mc.get()))], int(eval(Entrada_num_ensaios_mc.get())))
        grafico.plot(self.t, self.f(self.t,self.est_mc[0]), '-', label="Monte Carlo")
        grafico.legend()
        canvas.draw_idle()
        print("Kef_mc =", self.est_mc[0])
        Label_Kef_mc.grid(column=1, row=1, stick=tk.W, padx=5)
        Label_Kef_mc["text"]="Monte Carlo\nKef ={0} J/m³\nQuiquadrado ={1}\n ".format(self.est_mc[0], min(self.listaquis))
        janela.update()
        
    def ag_fit(self):
        self.est_ag, self.lista_ag, self.listaquis_ag = self.AG(self.f, self.xs, self.ys, [float(eval(Entrada_p0_ag.get()))], [float(eval(valor_raio_ag.get()))], int(eval(Entrada_limite.get())))
        grafico.plot(self.t, self.f(self.t,self.est_ag[0]), '-', label="Algoritmo Genético")
        grafico.legend()
        canvas.draw_idle()
        print("Kef_ag =", self.est_ag[0])
        Label_Kef_ag.grid(column=2, row=1, stick=tk.W, padx=5)
        Label_Kef_ag["text"]="Algoritmo Genético\nKef ={0} J/m³\nQuiquadrado ={1}\n ".format(self.est_ag[0], min(self.listaquis_ag))
        janela.update()
        
    def Ajuste(self):
        if var_curve_fit.get()==1:
            self.curve_fit()
        if var_MonteCarlo.get()==1:
            self.MonteCarlo_fit()
        if var_AG.get()==1:
            self.ag_fit()
        
    def mostrar_ag(self):
        if var_AG.get()==1:
            Label_p0_ag.grid(row=1, column=2)
            Entrada_p0_ag.grid(row=2, column=2)
            Label_raio_ag.grid(row=3, column=2)
            Entrada_raio_ag.grid(row=4, column=2)
            Label_limite.grid(row=5,column=2)
            Entrada_limite.grid(row=6, column=2)
        if var_AG.get()==0:
            Label_p0_ag.grid_forget()
            Entrada_p0_ag.grid_forget()
            Label_raio_ag.grid_forget()
            Entrada_raio_ag.grid_forget()
            Label_limite.grid_forget()
            Entrada_limite.grid_forget()
    
    def mostrar_mc(self):
        if var_MonteCarlo.get()==1:
            Label_p0_mc.grid(row=1, column=1)
            Entrada_p0_mc.grid(row=2, column=1)
            Label_raio_mc.grid(row=3, column=1)
            Entrada_raio_mc.grid(row=4, column=1)
            Label_num_ensaios_mc.grid(row=5,column=1)
            Entrada_num_ensaios_mc.grid(row=6, column=1)
        if var_MonteCarlo.get()==0:
            Label_p0_mc.grid_forget()
            Entrada_p0_mc.grid_forget()
            Label_raio_mc.grid_forget()
            Entrada_raio_mc.grid_forget()
            Label_num_ensaios_mc.grid_forget()
            Entrada_num_ensaios_mc.grid_forget()
            
    def mostrar_c(self):
        if var_curve_fit.get()==1:
            Label_p0.grid(row=1, column=0)
            Entrada_p0.grid(row=2, column=0)
        if var_curve_fit.get()==0:
            Label_p0.grid_forget()
            Entrada_p0.grid_forget()
            
    def Escrever(self):
        filename = tk.filedialog.askopenfilename()
        dados=pd.read_excel(filename)
        self.ts=dados["t"].values
        self.ms=dados["m"].values
        for i in range(0,len(self.ts)):
            dadosdig.insert(tk.END, str(self.ts[i])+'\t'+str(self.ms[i])+"\n")
        print("Os dados do arquivo inserido foram escritos")

c=A()

#Janela
janela=tk.Tk()
janela.title("MNC")
#janela.geometry(f'{1280}x{1024}')
janela.state('zoomed')
#janela.attributes('-fullscreen',True)
#w, h = janela.winfo_screenwidth(), janela.winfo_screenheight()
#janela.geometry("%dx%d+0+0" % (w, h))
janela["bg"]="grey94"

#Menus
barra_menu=tk.Menu(janela)
menu_Arquivo=tk.Menu(barra_menu, tearoff=0)
menu_Arquivo.add_command(label="Abrir", command=c.Escrever)
barra_menu.add_cascade(label="Arquivo", menu=menu_Arquivo)

janela.config(menu=barra_menu)

#Frames
frame_toolbar=tk.Frame(janela)
frame_toolbar.grid(row=3, column=2, columnspan=2, stick=tk.SW)

frame_plotagem=tk.LabelFrame(janela, text="Plotagem",bg="grey94")
frame_plotagem.grid(row=1, column=1, stick="n")

frame_param=tk.LabelFrame(janela, text="Parâmetros", bg="grey94")
frame_param.grid(row=0, column=1, stick="n")

frame_dados=tk.LabelFrame(janela, text="Dados", bg="grey94")
frame_dados.grid(row=0, column=0, rowspan=2, padx=5)

frame_saida=tk.LabelFrame(janela, text="Saída", bg="grey94", height=80)
frame_saida.grid(row=2, column=0, columnspan=2, rowspan=2, padx=5, stick="nsew")

#Figura onde será plotado o gráfico
fig=Figure(figsize=(6,6), dpi=100)
fig.suptitle("Massa aparente em função do tempo (normalizada)") 
    # fig.text(0.52,0.04, "Tempo (s)", ha="center", va="center")
    # fig.text(0.05,0.5, "Massa normalizada", ha="center", va="center", rotation=90)
fig.patch.set_facecolor('#F0F0F0')
grafico=fig.add_subplot(111)
grafico.set_xlabel("Tempo (s)")
grafico.set_ylabel("Massa normalizada")

canvas=FigureCanvasTkAgg(fig,master=janela)
canvas.draw()
canvas.get_tk_widget().grid(row=0, column=2, columnspan=2, rowspan=3, stick="nsew")

toolbar=NavigationToolbar2Tk(canvas,frame_toolbar)
toolbar.update()

#Variáveis
var_curve_fit=tk.IntVar()

var_MonteCarlo=tk.IntVar()

var_AG=tk.IntVar()

var_barra=tk.IntVar()

valor_mt=tk.StringVar()

valor_rho=tk.StringVar()

valor_uim=tk.StringVar()
valor_uim.set("0.453")

valor_g=tk.StringVar()
valor_g.set("9.82")

valor_z0=tk.StringVar()
valor_z0.set("0.0063")

valor_z1=tk.StringVar()

valor_R=tk.StringVar()
valor_R.set("0.0049")

valor_Msv=tk.StringVar()

valor_Dm=tk.StringVar()

valor_sigmad=tk.StringVar()

valor_u0=tk.StringVar()
valor_u0.set("1.256*10**-6")

valor_H=tk.StringVar()

valor_Kb=tk.StringVar()
valor_Kb.set("1.3807*10**-23")

valor_T=tk.StringVar()
valor_T.set("300")

valor_t0=tk.StringVar()
valor_t0.set("10**-9")

valor_p0=tk.StringVar()

valor_p0_mc=tk.StringVar()

valor_p0_ag=tk.StringVar()

valor_raio_mc=tk.StringVar()

valor_raio_ag=tk.StringVar()

valor_num_ensaios=tk.StringVar()
valor_num_ensaios.set('100')

valor_limite=tk.StringVar()

#Caixas de texto
Entrada_mt=tk.Entry(frame_param, width=21, textvariable=valor_mt)
Entrada_mt.grid(row=2, column=1,pady=3,padx=13)

Entrada_rho=tk.Entry(frame_param,width=21, textvariable=valor_rho)
Entrada_rho.grid(row=3, column=1,pady=3,padx=13)

Entrada_uim=tk.Entry(frame_param,width=21, textvariable=valor_uim)
Entrada_uim.grid(row=4, column=1,pady=3,padx=13)

Entrada_g=tk.Entry(frame_param,width=21, textvariable=valor_g)
Entrada_g.grid(row=5, column=1,pady=3,padx=13)

Entrada_z0=tk.Entry(frame_param,width=21, textvariable=valor_z0)
Entrada_z0.grid(row=6, column=1,pady=3,padx=13)

Entrada_z1=tk.Entry(frame_param,width=21, textvariable=valor_z1)
Entrada_z1.grid(row=7, column=1,pady=3,padx=13)

Entrada_R=tk.Entry(frame_param,width=21, textvariable=valor_R)
Entrada_R.grid(row=8, column=1,pady=3,padx=13)

Entrada_Msv=tk.Entry(frame_param,width=21, textvariable=valor_Msv)
Entrada_Msv.grid(row=9, column=1,pady=3,padx=13)

Entrada_Dm=tk.Entry(frame_param,width=21, textvariable=valor_Dm)
Entrada_Dm.grid(row=10, column=1,pady=3,padx=13)

Entrada_sigmad=tk.Entry(frame_param,width=21, textvariable=valor_sigmad)
Entrada_sigmad.grid(row=11, column=1,pady=3,padx=13)

Entrada_u0=tk.Entry(frame_param,width=21, textvariable=valor_u0)
Entrada_u0.grid(row=12, column=1,pady=3,padx=13)

Entrada_H=tk.Entry(frame_param,width=21, textvariable=valor_H)
Entrada_H.grid(row=13, column=1,pady=3,padx=13)

Entrada_Kb=tk.Entry(frame_param,width=21, textvariable=valor_Kb)
Entrada_Kb.grid(row=14, column=1,pady=3,padx=13)

Entrada_T=tk.Entry(frame_param,width=21, textvariable=valor_T)
Entrada_T.grid(row=15, column=1,pady=3,padx=13)

Entrada_t0=tk.Entry(frame_param,width=21, textvariable=valor_t0)
Entrada_t0.grid(row=16, column=1,pady=3,padx=13)

Entrada_p0=tk.Entry(frame_plotagem,width=21, textvariable=valor_p0)

Entrada_p0_mc=tk.Entry(frame_plotagem,width=21, textvariable=valor_p0_mc)

Entrada_p0_ag=tk.Entry(frame_plotagem,width=21, textvariable=valor_p0_ag)

Entrada_num_ensaios_mc=tk.Entry(frame_plotagem, width=21, textvariable=valor_num_ensaios)

Entrada_raio_mc=tk.Entry(frame_plotagem,width=21, textvariable=valor_raio_mc)

Entrada_raio_ag=tk.Entry(frame_plotagem,width=21, textvariable=valor_raio_ag)

Entrada_limite=tk.Entry(frame_plotagem,width=21, textvariable=valor_limite)

#Labels
Label_dretorio=tk.Label(frame_dados, width=40, text="Escreva os dados obtidos em duas colunas separadas\npor tab: (uma coluna para tempo e outra para massa)", fg="black")
Label_dretorio.grid(row=0, column=0, columnspan=2)

Label_mt=tk.Label(frame_param, text="Massa total (kg)", fg="black")
Label_mt.grid(row=2, column=0, padx=13)

Label_rho=tk.Label(frame_param, text="Densidade da amostra (kg/m)",fg="black")
Label_rho.grid(row=3, column=0, padx=13)

Label_uim=tk.Label(frame_param, text="Momento Magnético do Ímã (A.m²)",fg="black")
Label_uim.grid(row=4, column=0, padx=13)

Label_g=tk.Label(frame_param, text="Aceleração da gravidade (m/s²)",fg="black")
Label_g.grid(row=5, column=0, padx=13)

Label_z0=tk.Label(frame_param, text="z0 (m)",fg="black")
Label_z0.grid(row=6, column=0, padx=13)

Label_z1=tk.Label(frame_param, text="z1 (m)",fg="black")
Label_z1.grid(row=7, column=0, padx=13)

Label_R=tk.Label(frame_param, text="Raio do porta amostra (m)",fg="black")
Label_R.grid(row=8, column=0, padx=13)

Label_Msv=tk.Label(frame_param, text="Magnetização de Saturação Volumétrica (A/m)",fg="black")
Label_Msv.grid(row=9, column=0, padx=13)

Label_Dm=tk.Label(frame_param, text="Diâmetro médio das nanopartículas (m)",fg="black")
Label_Dm.grid(row=10, column=0, padx=13)

Label_sigmad=tk.Label(frame_param, text="Fator de dispersão do Diâmetro",fg="black")
Label_sigmad.grid(row=11, column=0, padx=13)

Label_u0=tk.Label(frame_param, text="Permeabilidade Magnética do meio (N/A²)",fg="black")
Label_u0.grid(row=12, column=0, padx=13)

Label_H=tk.Label(frame_param, text="Campo Magnético (A/m)",fg="black")
Label_H.grid(row=13, column=0, padx=13)

Label_Kb=tk.Label(frame_param, text="Constante de Boltzman (m².kg/K.s)",fg="black")
Label_Kb.grid(row=14, column=0, padx=13)

Label_T=tk.Label(frame_param, text="Temperatura (K)",fg="black")
Label_T.grid(row=15, column=0, padx=13)

Label_t0=tk.Label(frame_param, text="Fator de frequência (s)",fg="black")
Label_t0.grid(row=16, column=0, padx=13)

Label_p0=tk.Label(frame_plotagem, text="Estimativa de Kef (J/m³)", fg="black")

Label_p0_mc=tk.Label(frame_plotagem, text="Estimativa de Kef (J/m³)", fg="black")

Label_p0_ag=tk.Label(frame_plotagem, text="Estimativa de Kef (J/m³)", fg="black")

Label_raio_mc=tk.Label(frame_plotagem, text="Raio", fg="black")

Label_raio_ag=tk.Label(frame_plotagem, text="Raio", fg="black")

Label_limite=tk.Label(frame_plotagem, text="Limite", fg="black")

Label_num_ensaios_mc=tk.Label(frame_plotagem, text="Nº de ensaios", fg="black")

Label_Kef_c=tk.Label(frame_saida, text="")
Label_Kef_mc=tk.Label(frame_saida, text="")
Label_Kef_ag=tk.Label(frame_saida, text="")

#Label_Ensaios=tk.Label(frame_saida)

#Text
scrollbar1 = tk.Scrollbar(frame_dados)
scrollbar1.grid(row=3, column=3, stick="ns")
dadosdig=tk.Text(frame_dados, width=35, height=33, yscrollcommand=scrollbar1.set)
dadosdig.grid(row=3, column=0, columnspan=2, pady=5)
scrollbar1.config(command = dadosdig.yview)

#CheckButtons
check_curve_fit=tk.Checkbutton(frame_plotagem, text="Min. Quadrados", variable=var_curve_fit, command=c.mostrar_c)
check_curve_fit.grid(row=0, column=0, padx=15)

check_MonteCarlo=tk.Checkbutton(frame_plotagem, text="Monte Carlo", variable=var_MonteCarlo, command=c.mostrar_mc)
check_MonteCarlo.grid(row=0, column=1, padx=15)

check_AG=tk.Checkbutton(frame_plotagem, text="Algoritmo Genético", variable=var_AG, command=c.mostrar_ag)
check_AG.grid(row=0, column=2, padx=15)

#Botões
BotPlot=tk.Button(frame_plotagem,text="Plotar",command=c.Variaveis,fg="black", width=40, height=1).grid(row=7, column=0, columnspan=3, pady=5, padx=5)

janela.mainloop()