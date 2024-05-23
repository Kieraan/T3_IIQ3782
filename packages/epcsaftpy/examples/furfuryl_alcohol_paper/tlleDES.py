import numpy as np
from sgtpy.equilibrium import lle, lle_init          # sgtpy module
import pandas as pd
import matplotlib.pyplot as plt
import ternary
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell

def mol_to_mass(x, saft):
    return x * saft.Mw / np.dot(x, saft.Mw)

def mass_to_mol(w, saft):
    return (w / saft.Mw) / np.dot(w, 1/saft.Mw)

def solving_ternary_DES(T, P, s, x0, w0, z0, zf, n, saft, new_mode):
    dz = (zf - z0)/n
    
    # List to store the information
    x, w, z = [], [], []
    
    # Solving the coexistence curve
    v0 = [None, None]
    Xass0 = [None, None]
    l = 0
    for i in range(2*n):
        out = lle(x0, w0, z0, T, P, saft, 
                  v0 = v0, Xass0 = Xass0, full_output = True, nacc = 10)
        if abs(out.error_inner)<1e-5 and abs(out.error_outer<1e-5) and abs(out.v[0] - out.v[1])>1e-8:
            x0, w0 = out.X
            v0 = out.v
            Xass0 = out.Xass
            z.append(mol_to_mass(z0, saft)), x.append(mol_to_mass(x0, saft)), w.append(mol_to_mass(w0, saft))
            if new_mode:
                z0 = (x0 + w0)/2
                z_DES = z0[2] + z0[3]
                z0 = np.array([z0[0], z0[1], z_DES * s/(s + 1), z_DES / (s + 1)])
            z0 += dz
            while min(z0)<0:
                z0 -= dz
                dz /= 2.0
                z0 += dz
                l += 1
                if l > 10:
                    break
        else:
            z0 -= dz
            dz /= 2.0
            z0 += dz
            l += 1
            if l > 10:
                break
    return z, x, w


def lle_DES(saft, T, P, s, name_DES, new_mode = False, n = 200, zW0 = 0.8, zWf = 0.85):
    data = pd.read_excel('data_exp/LLE_water_FA_DES.xlsx', name_DES)
    w1E, w2E, w3E = np.array(data.iloc[:,[0,1,2]]).T
    x1E, x2E, x3E = np.array(data.iloc[:,[3,4,5]]).T

    # Initial point
    zFA = 0.0
    zDES = 1 - zW0
    zHBD = s * zDES/(s + 1)
    zHBA = zDES/(s + 1)
    z0 = np.array([zW0, zFA, zHBD, zHBA])
    x0, w0 = lle_init(z0, T, P, saft)
    
    # Final point
    zFA = 1 - zWf
    zf = np.array([zWf, zFA, 0., 0.])
    
    # Solving coexistence
    z, x, w = solving_ternary_DES(T, P, s, x0, w0, z0, zf, n, saft, new_mode)          
                
    # Solving for the experimental tie-line        
    nE = len(w1E)
    xT, wT = [], []
    for i in range(nE):
        x0 = mass_to_mol(np.array([x1E[i], x2E[i], x3E[i] * s / (s + 1), x3E[i]  / (s + 1)]), saft)
        w0 = mass_to_mol(np.array([w1E[i], w2E[i], w3E[i] * s / (s + 1), w3E[i]  / (s + 1)]), saft)
        z0 = (x0 + w0)/2
        z_DES = z0[2] + z0[3]
        z0 = np.array([z0[0], z0[1], z_DES * s/(s + 1), z_DES / (s + 1)])
        x0, w0 = lle_init(z0, T, P, saft)
        out = lle(x0, w0, z0, T, P, saft, full_output = True)
        if abs(out.error_inner)<1e-5 and abs(out.error_outer<1e-5) and abs(out.v[0] - out.v[1])>1e-8:
            x0, w0 = out.X
            xT.append(mol_to_mass(x0, saft)), wT.append(mol_to_mass(w0, saft))
    
    return z, x, w, xT, wT


def show_plot(z, x, w, xT, wT, name_DES, name_option, show_z = False):
    data = pd.read_excel('data_exp/LLE_water_FA_DES.xlsx', name_DES)
    w1E, w2E, w3E = np.array(data.iloc[:,[0,1,2]]).T
    x1E, x2E, x3E = np.array(data.iloc[:,[3,4,5]]).T
    data = pd.read_excel('data_exp/LLE_water_FA_DES.xlsx', name_DES + "_SIM")
    w1S, w2S, w3S = np.array(data.iloc[:,[0,1,2]]).T
    x1S, x2S, x3S = np.array(data.iloc[:,[3,4,5]]).T
    
    
    colorSAFT = "#49C60A"
    colorEXP = "#FC7725"
    colorSIM = "#404080"
    
    font = {'weight' : 'normal',
            'size'   : 16}
    plt.rc('font', **font)
    plt.rcParams['figure.dpi'] = 300
        
        
    figure, tax = ternary.figure(scale = 1.0)
    figure.set_size_inches(8, 8)
    
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="black", multiple=0.1, linewidth=0.5)
    tax.ticks(axis='lbr', linewidth=1, multiple=1)
    
    # Simulation tie-line
    nS = len(x1S)
    for i in range(nS):
        p1 = (x1S[i], x2S[i], x3S[i])
        p2 = (w1S[i], w2S[i], w3S[i])
        tax.line(p1, p2, linewidth = 1.5, markersize = 8, 
                 markeredgewidth = 1, markeredgecolor = "black", 
                 marker = "X", color = colorSIM, linestyle = ":")
        
    
    # Experimental tie-line
    nE = len(x1E)
    for i in range(nE):
        p1 = (x1E[i], x2E[i], x3E[i])
        p2 = (w1E[i], w2E[i], w3E[i])
        tax.line(p1, p2, linewidth = 1.5, markersize = 8, 
                 markeredgewidth = 1, markeredgecolor = "black", 
                 marker = "o", color = colorEXP, linestyle = ":")
        
    # PC-SAFT tie-line
    nT = len(xT)
    for i in range(nT):
        p1 = (xT[i][0], xT[i][1], xT[i][2] + xT[i][3])
        p2 = (wT[i][0], wT[i][1], wT[i][2] + wT[i][3])
        tax.line(p1, p2, linewidth = 1.5, markersize = 8, 
                 markeredgewidth = 1, markeredgecolor = "black", 
                 marker = "s", color = colorSAFT, linestyle = ":")
        
    # PC-SAFT plot
    n = len(x)
    X, W, Z = np.zeros([n, 3]), np.zeros([n, 3]), np.zeros([n, 3])
    for i in range(n):
        X[i, :] = x[i][0], x[i][1], x[i][2] + x[i][3]
        W[i, :] = w[i][0], w[i][1], w[i][2] + w[i][3]
        Z[i, :] = z[i][0], z[i][1], z[i][2] + z[i][3]
    tax.plot(X, color = colorSAFT)
    tax.plot(W, color = colorSAFT)
    if show_z:
        tax.plot(Z, color = colorSAFT)
    
    # Plot options
    tax.ticks(axis='lbr', multiple=0.1, linewidth=1, tick_formats='%.1f')
    fontsize = 15
    tax.left_axis_label( name_DES + " (3)", fontsize=fontsize, offset = 0.1)
    tax.right_axis_label('Furfuryl Alcohol (2)', fontsize=fontsize, offset = 0.1)
    tax.bottom_axis_label('Water (1)', fontsize=fontsize, offset = -0.05)
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    
    plt.axis('off')
    plt.savefig('figures/LLE_W+FA+'+ name_DES + ' (' + name_option + ').pdf')
    plt.show()
        

    
    Excel = xlsxwriter.Workbook('results/LLE_W+FA+'+ name_DES 
                                + ' (' + name_option +  ').xlsx')
    Hoja1 = Excel.add_worksheet("LLE")
    l = 0
    Hoja1.write(xl_rowcol_to_cell(l, 0),"Exp - Tie line")
    Hoja1.write(xl_rowcol_to_cell(l, 3),"MD - Tie line")
    Hoja1.write(xl_rowcol_to_cell(l, 6),"PC-SAFT - Tie line")
    Hoja1.write(xl_rowcol_to_cell(l, 9),"PC-SAFT - DES-phase")
    Hoja1.write(xl_rowcol_to_cell(l, 12),"PC-SAFT - Water-phaser")
    l += 1
    Hoja1.write(xl_rowcol_to_cell(l, 0),"Exp-Water")
    Hoja1.write(xl_rowcol_to_cell(l, 1),"Exp-FA")
    Hoja1.write(xl_rowcol_to_cell(l, 2),"Exp-DES")
    Hoja1.write(xl_rowcol_to_cell(l, 3),"MD-Water")
    Hoja1.write(xl_rowcol_to_cell(l, 4),"MD-FA")
    Hoja1.write(xl_rowcol_to_cell(l, 5),"MD-DES")
    Hoja1.write(xl_rowcol_to_cell(l, 6),"SAFT-TL-Water")
    Hoja1.write(xl_rowcol_to_cell(l, 7),"SAFT-TL-FA")
    Hoja1.write(xl_rowcol_to_cell(l, 8),"SAFT-TL-DES")
    Hoja1.write(xl_rowcol_to_cell(l, 9),"SAFT-Water")
    Hoja1.write(xl_rowcol_to_cell(l, 10),"SAFT-FA")
    Hoja1.write(xl_rowcol_to_cell(l, 11),"SAFT-DES")    
    Hoja1.write(xl_rowcol_to_cell(l, 12),"SAFT-Water")
    Hoja1.write(xl_rowcol_to_cell(l, 13),"SAFT-FA")
    Hoja1.write(xl_rowcol_to_cell(l, 14),"SAFT-DES")
    # Exp - Tie line
    n = len(w1E)
    l = 1
    for i in range(n):
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 0), w1E[i])
        Hoja1.write(xl_rowcol_to_cell(l, 1), w2E[i])
        Hoja1.write(xl_rowcol_to_cell(l, 2), w3E[i])
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 0), x1E[i])
        Hoja1.write(xl_rowcol_to_cell(l, 1), x2E[i])
        Hoja1.write(xl_rowcol_to_cell(l, 2), x3E[i])
        l += 1
    # MD - Tie line
    n = len(w1S)
    l = 1
    for i in range(n):
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 3), w1S[i])
        Hoja1.write(xl_rowcol_to_cell(l, 4), w2S[i])
        Hoja1.write(xl_rowcol_to_cell(l, 5), w3S[i])
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 3), x1S[i])
        Hoja1.write(xl_rowcol_to_cell(l, 4), x2S[i])
        Hoja1.write(xl_rowcol_to_cell(l, 5), x3S[i])
        l += 1
    # SAFT - Tie line
    n = len(xT[:][0])
    l = 1
    for i in range(n):
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 6), wT[i][0])
        Hoja1.write(xl_rowcol_to_cell(l, 7), wT[i][1])
        Hoja1.write(xl_rowcol_to_cell(l, 8), wT[i][2] + wT[i][3])
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 6), xT[i][0])
        Hoja1.write(xl_rowcol_to_cell(l, 7), xT[i][1])
        Hoja1.write(xl_rowcol_to_cell(l, 8), xT[i][2] + xT[i][3])
        l += 1

    # LLE
    n = len(W[:,0])
    l = 1
    for i in range(n):
        l += 1
        Hoja1.write(xl_rowcol_to_cell(l, 9), W[i,0])
        Hoja1.write(xl_rowcol_to_cell(l, 10), W[i,1])
        Hoja1.write(xl_rowcol_to_cell(l, 11), W[i,2])
        Hoja1.write(xl_rowcol_to_cell(l, 12), X[i,0])
        Hoja1.write(xl_rowcol_to_cell(l, 13), X[i,1])
        Hoja1.write(xl_rowcol_to_cell(l, 14), X[i,2])
        
        
    Hoja2 = Excel.add_worksheet("D")
    Hoja3 = Excel.add_worksheet("S")
    l = 0
    Hoja2.write(xl_rowcol_to_cell(l, 0),"Exp")
    Hoja2.write(xl_rowcol_to_cell(l, 2),"MD")
    Hoja2.write(xl_rowcol_to_cell(l, 4),"PC-SAFT")
    
    Hoja3.write(xl_rowcol_to_cell(l, 0),"Exp")
    Hoja3.write(xl_rowcol_to_cell(l, 2),"MD")
    Hoja3.write(xl_rowcol_to_cell(l, 4),"PC-SAFT")
    l += 1
    Hoja2.write(xl_rowcol_to_cell(l, 0),"xFA in W")
    Hoja2.write(xl_rowcol_to_cell(l, 1),"D")
    Hoja2.write(xl_rowcol_to_cell(l, 2),"xFA in W")
    Hoja2.write(xl_rowcol_to_cell(l, 3),"D")    
    Hoja2.write(xl_rowcol_to_cell(l, 4),"xFA in W")
    Hoja2.write(xl_rowcol_to_cell(l, 5),"D")
    
    Hoja3.write(xl_rowcol_to_cell(l, 0),"xFA in W")
    Hoja3.write(xl_rowcol_to_cell(l, 1),"S")
    Hoja3.write(xl_rowcol_to_cell(l, 2),"xFA in W")
    Hoja3.write(xl_rowcol_to_cell(l, 3),"S")    
    Hoja3.write(xl_rowcol_to_cell(l, 4),"xFA in W")
    Hoja3.write(xl_rowcol_to_cell(l, 5),"S")
    # Exp
    n = len(w1E)
    l = 1
    for i in range(n):
        if x2E[i] != 0:
            l += 1
            with np.errstate(divide='ignore', invalid='ignore'):
                Dexp = w2E[i] / x2E[i]
                Sexp = Dexp / (w1E[i] / x1E[i])
            Hoja2.write(xl_rowcol_to_cell(l, 0), x2E[i])
            Hoja2.write(xl_rowcol_to_cell(l, 1), Dexp)
            Hoja3.write(xl_rowcol_to_cell(l, 0), x2E[i])
            Hoja3.write(xl_rowcol_to_cell(l, 1), Sexp)
    # MD
    n = len(w1S)
    l = 1
    for i in range(n):
        if w2S[i] != 0:
            l += 1
            with np.errstate(divide='ignore', invalid='ignore'):
                DS = w2S[i] / x2S[i]
                SS = DS / (w1S[i] / x1S[i])
            Hoja2.write(xl_rowcol_to_cell(l, 2), x2S[i])
            Hoja2.write(xl_rowcol_to_cell(l, 3), DS)
            Hoja3.write(xl_rowcol_to_cell(l, 2), x2S[i])
            Hoja3.write(xl_rowcol_to_cell(l, 3), SS)
    # SAFT
    x1m, x2m =  np.array(x)[:,0], np.array(x)[:,1]
    w1m, w2m =  np.array(w)[:,0], np.array(w)[:,1]
    n = len(x1m)
    l = 1
    for i in range(n):
        if x2m[i] != 0:
            l += 1
            with np.errstate(divide='ignore', invalid='ignore'):
                Dsaft = w2m[i] / x2m[i]
                Ssaft = Dsaft / (w1m[i] / x1m[i])
            Hoja2.write(xl_rowcol_to_cell(l, 4), x2m[i])
            Hoja2.write(xl_rowcol_to_cell(l, 5), Dsaft)
            Hoja3.write(xl_rowcol_to_cell(l, 4), x2m[i])
            Hoja3.write(xl_rowcol_to_cell(l, 5), Ssaft)    
        
    Excel.close()


def insert_plot(name_DES, name_option, color, x, w, axs):
    data = pd.read_excel('data_exp/LLE_water_FA_DES.xlsx', name_DES)
    w1E, w2E, w3E = np.array(data.iloc[:,[0,1,2]]).T
    x1E, x2E, x3E = np.array(data.iloc[:,[3,4,5]]).T
    data = pd.read_excel('data_exp/LLE_water_FA_DES.xlsx', name_DES + "_SIM")
    w1S, w2S, w3S = np.array(data.iloc[:,[0,1,2]]).T
    x1S, x2S, x3S = np.array(data.iloc[:,[3,4,5]]).T
    
    x1m, x2m =  np.array(x)[:,0], np.array(x)[:,1]
    w1m, w2m =  np.array(w)[:,0], np.array(w)[:,1]

    with np.errstate(divide='ignore', invalid='ignore'):
        Dsaft = w2m/x2m
        Ssaft = (w2m/x2m)/(w1m/x1m)
        
        Dexp = w2E/x2E
        Sexp = (w2E/x2E)/(w1E/x1E)
        
        DS = w2S/x2S
        SS = (w2S/x2S)/(w1S/x1S)
    
    axs[0].plot(x2m, Ssaft, color)
    axs[1].plot(x2m, Dsaft, color, label = name_DES)
    axs[0].plot(x2E, Sexp, markersize = 10, markeredgewidth=1, markeredgecolor="black", c=color, marker="o", 
                linewidth=0)
    axs[1].plot(x2E, Dexp, markersize = 10, markeredgewidth=1, markeredgecolor="black", c=color, marker="o", 
                linewidth=0)
    axs[0].plot(x2S, SS, markersize = 10, markeredgewidth=1, markeredgecolor="black", c=color, marker="X", 
                linewidth=0)
    axs[1].plot(x2S, DS, markersize = 10, markeredgewidth=1, markeredgecolor="black", c=color, marker="X", 
                linewidth=0)
