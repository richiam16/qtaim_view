import sympy as sp   # For primitive expression
# import itertools as it  #Para realizar iteracciones de manera sencilla
from sympy.abc import x, y, z
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Surface Plots
from matplotlib import cm


class PG:
    """Class to definie primitive functions to be used in the construction of the density"""

    def __init__(self, centre, type, exp):
        """Centre: string, the # of centre within the .wfn file (# nuclei)
           type: str, the type of function to be used
           exp: float, the exponential of the gaussian part"""
        self.centre = float(centre)
        self.type = type
        self.exp = exp

    def cons_exp(self, CI):
        exout1 = "Obtaining Primitive Function from center" + \
            " "+str(self.centre)+" "+"type:"+" "+str(self.type)
        print(exout1)
        # myfile=open("/home/ricardo/Escritorio/Corr/QTAIM/typeG.txt","r") #abre archivo que contine informacion de A
        myfile = open("/mnt/c/Users/richi/Desktop/QTAIM/typeG.txt")
        lines = myfile.readlines()  # archiva las lineas del archivo
        x, y, z = sp.symbols("x y z")
        # cambia la D por E para exponentes
        self.exp = self.exp[:-4]+"E"+self.exp[len(self.exp)-3:]
        # convierte el exponente de un string a numero
        self.exp = float(self.exp)
        # El loop busca en el archivo el tipo de funcion y asigna exponentes a A
        for i in range(0, len(lines)):
            prueb = lines[i].split()
            if prueb[0] == self.type:
                fx = x**float(prueb[1])
                fy = y**float(prueb[2])
                fz = z**float(prueb[3])
                coefx, coefy, coefz = float(prueb[1]), float(
                    prueb[2]), float(prueb[3])
        self.A = fx*fy*fz  # Define la funcion A
        self.td = (((8*self.exp)**(coefx+coefy+coefz)*sp.factorial(coefx)*sp.factorial(coefy)*sp.factorial(coefz)) /
                   (sp.factorial(2*coefx)*sp.factorial(2*coefy)*sp.factorial(2*coefz)))**(1./2.)  # Para la funcion factorial
        # el factor que depende de pi
        self.ae = sp.N(((2*self.exp)/sp.pi)**(3./4.))
        for i in range(0, len(CI)):  # este loop saca las coordenadas del centro correspondiente
            if self.centre == CI[i][1]:
                self.cx, self.cy, self.cz = CI[i][2], CI[i][3], CI[i][4]
        # EScribe la funcion exponencial
        self.F = self.td*self.ae * \
            sp.exp(-self.exp*((x+self.cx)**2+(y+self.cy)**2+(z+self.cz)**2))
        self.G = self.F*self.A  # Escribre la primitiva total

    def primeG(self):
        """In construction, this function is intented to copy PrimeG function of Baders original code"""
        self.dev_x = sp.diff(self.G, "x")
        self.dev_y = sp.diff(self.G, "y")
        self.dev_z = sp.diff(self.G, "z")
        # self.dev_4=list(it.combinations_with_replacement([x,y,z],4))
        # self.dev_3=[self.dev_4[i][:-1] for i in range(0,len(self.dev_4))]
        # self.dev_2=[self.dev_3[i][:-1] for i in range(0,len(self.dev_3))]
        # self.dev_1=[self.dev_2[i][:-1] for i in range(0,len(self.dev_2))]
        # self.DGG=[]
        # for i in range(0,len(self.dev_4)):
        # dev_par_1=sp.diff(self.F,self.dev_4[i][0])
        # dev_par_2=sp.diff(dev_par_1,self.dev_4[i][1])
        # dev_par_3=sp.diff(dev_par_2,self.dev_4[i][2])
        # dev_par_4=sp.diff(dev_par_3,self.dev_4[i][3])
        # self.DGG.append([dev_par_1,dev_par_2,dev_par_3,dev_par_4])

        # self.dG4=[self.DGG[i][3:] for i in range(0,len(self.DGG))]
        # self.dG3=[self.DGG[0][2],self.DGG[3][2],self.DGG[5][2],self.DGG[6][2],self.DGG[8][2],self.DGG[9][2],self.DGG[10][2],self.DGG[12][2],self.DGG[13][2],self.DGG[14][2]]
        # self.dG2=[self.DGG[0][1],self.DGG[6][1],self.DGG[9][1],self.DGG[10][1],self.DGG[13][1],self.DGG[14][1]]
        # self.dG1=[self.DGG[0][0],self.DGG[10][0],self.DGG[14][0]]


class Mol:
    """Class to define the density of molecule"""

    def __init__(self, path):
        """path:  str, path to find the .wfn file"""
        self.path = path

    def rdwfn(self):
        """Method of read the contents of the .wfn file"""
        print("Reading File")
        myfile = open(self.path, "r")
        lines = myfile.readlines()
        self.nom = lines[0][:-1]  # Takes the filename of the .wfn file
        temp = lines[1].split()
        self.btype = temp[0]  # Type of primitive functions
        self.num_mol_orn = float(temp[1])  # Total number of molecular orbitals
        # Number number of primitive gaussian functions
        self.num_prim = float(temp[4])
        self.num_nucl = float(temp[6])  # Number of nuclei in the file
        self.CI = []  # Coordinate Information (xyz)
        self.M = []  # Molecular orbitals
        self.MC = []  # Coefficient of the primitive functions
        centre = []
        type = []
        exp = []
        # de la linea dos hasta el numero de nucleos
        for i in range(2, 2+int(self.num_nucl)):
            temp_2 = lines[i].split()
            self.CI.append([temp_2[0], float(temp_2[1]), float(temp_2[4]), float(temp_2[5]), float(
                temp_2[6]), float(temp_2[9])])  # Informacion: Nucleo, #, x,y,z,carga
        # del final de las lineas de nucleos, hasta el fin del documento
        for i in range(2+int(self.num_nucl), len(lines)):
            clasi = lines[i].split()
            if clasi[0] == "CENTRE":
                centre = centre+clasi[2:]
            if clasi[0] == "TYPE":
                type = type+clasi[2:]
            if clasi[0] == "EXPONENTS":
                exp = exp+clasi[1:]
            if clasi[0] == "MO":
                self.M.append([clasi[1], clasi[7], clasi[11]])
                # self.M.append([clasi[1],clasi[5],clasi[8]])
                conf = 0
                OM = []
                for j in range(i, len(lines)):
                    if len(OM) != self.num_prim:
                        conf = conf+1
                        app = lines[conf+i].split()
                        OM = OM+app
                self.MC.append(OM)

            if clasi[0] == "END":
                temp_4 = lines[i+1].split()
                self.ET = float(temp_4[3])
                # self.ET=float(temp_4[4])
                self.viral = float(temp_4[6])
                # self.viral=float(temp_4[7])

        self.GF = [PG(centre[i], type[i], exp[i])
                   for i in range(0, len(centre))]
        print("End of file")

    def gaus4(self):
        """Method to reconstruct density from gaussian primitives"""
        self.rho = []
        for i in range(0, len(self.GF)):
            self.GF[i].cons_exp(self.CI)
            self.GF[i].primeG()
        print("Obtaining MOs")
        for i in range(0, len(self.MC)):
            psi = 0
            for j in range(0, len(self.MC[0])):
                self.MC[i][j] = self.MC[i][j][:-4]+"E" + \
                    self.MC[i][j][len(self.MC[i][j])-3:]
                self.MC[i][j] = float(self.MC[i][j])
                psi = psi+self.MC[i][j]*self.GF[j].G
            self.rho.append([psi*psi])
            # An f string would work better here"
            print("Done MO"+" "+str(i+1)+" "+"of"+" "+str(self.num_mol_orn))
        print("Rho is ready")
        self.rho_an = self.rho[0]
        # for i in range(0,len(self.rho)):
        #        self.rho[i]=sp.lambdify([x,y,z],self.rho[i],"numpy")
        # print("Rho is now a function")

    def graph(self, en=0, xyz="xy", min=0, max=0, pun=100, lim_rho=1, cont=0, min_c=30):
        """Method to graph the molecular orbitals densities, both isodensities and surface
        en: int, Number of molecular orbital, from 0 to ...
        XYZ: str, position of the graph, posible values: "x","y","z","xy","xz","yz"
        min: float, minimum xyz value of where to start the graph (in Angstroms)
        max: float, maximun xyz value of where to end the graph (in Angstroms)
        pun: float, number of points to evaluate in self.rho[en]
        lim_rho: float, limit value of rho allowed (max values=1)
        cont: int, isodensities or surface; 1 isodensities, anyother number surfaces
        min_c: number of isodensities"""
        rho_prub = sp.lambdify([x, y, z], self.rho[int(en)], "numpy")
        val_x = np.linspace(float(min), float(max), int(pun))
        val_y = np.linspace(float(min), float(max), int(pun))
        val_z = np.linspace(float(min), float(max), int(pun))
        fig = plt.figure()
        if xyz == "z":
            grap = rho_prub(0, 0, val_z)
            plt.plot(val_z, grap[0])
        if xyz == "y":
            grap = rho_prub(0, val_y, 0)
            plt.plot(val_y, grap[0])
        if xyz == "x":
            grap = rho_prub(val_x, 0, 0)
            plt.plot(val_x, grap[0])
        if xyz == "xy":
            val_x, val_y = np.meshgrid(val_x, val_y)
            self.grap = rho_prub(val_x, val_y, 0)
            for i in range(0, len(self.grap[0])):
                for j in range(0, len(self.grap[0][i])):
                    if self.grap[0][i][j] > float(lim_rho):
                        self.grap[0][i][j] = float(lim_rho)
            if float(cont) == 1:
                for i in range(0, len(self.grap[0])):
                    for j in range(0, len(self.grap[0][i])):
                        self.grap[0][i][j] = np.log10(self.grap[0][i][j])

                ax = fig.add_subplot(111)
                cs = ax.contour(val_x, val_y, self.grap[0], int(
                    min_c), cmap=cm.coolwarm)
            else:
                # ax=fig.gca(projection="3d")
                ax = fig.add_subplot(projection='3d')
                surf = ax.plot_surface(
                    val_x, val_y, self.grap[0], cmap=cm.coolwarm)
        if xyz == "xz":
            val_x, val_z = np.meshgrid(val_x, val_z)
            grap = rho_prub(val_x, 0, val_z)
            for i in range(0, len(grap[0])):
                for j in range(0, len(grap[0][i])):
                    if grap[0][i][j] > float(lim_rho):
                        grap[0][i][j] = float(lim_rho)
            if float(cont) == 1:
                for i in range(0, len(grap[0])):
                    for j in range(0, len(grap[0][i])):
                        grap[0][i][j] = np.log10(grap[0][i][j])

                ax = fig.add_subplot(111)
                cs = ax.contour(val_x, val_z, grap[0], int(min_c))
            else:
                #ax = fig.gca(projection="3d")
                ax = fig.add_subplot(projection='3d')
                #ax_1 = fig.gca(projection="3d")
                ax = fig.add_subplot(projection='3d')
                grad = np.gradient(grap[0])
                surf = ax.plot_surface(val_x, val_z, grap[0], cmap=cm.coolwarm)

        if xyz == "yz":
            # ax = fig.gca(projection="3d")
            ax = fig.add_subplot(projection='3d')
            val_y, val_z = np.meshgrid(val_y, val_z)
            grap = rho_prub(0, val_y, val_z)
            for i in range(0, len(grap[0])):
                for j in range(0, len(grap[0][i])):
                    if grap[0][i][j] > float(lim_rho):
                        grap[0][i][j] = float(lim_rho)
            if float(cont) == 1:
                for i in range(0, len(grap[0])):
                    for j in range(0, len(grap[0][i])):
                        grap[0][i][j] = np.log10(grap[0][i][j])

                ax = fig.add_subplot(111)
                cs = ax.contour(val_y, val_z, grap[0], int(min_c))
            else:
                surf = ax.plot_surface(val_y, val_z, grap[0], cmap=cm.coolwarm)
        # plt.axis('off')
        plt.show()

def main(path):
        density = Mol(path)
        density.rdwfn()
        density.gaus4()
        return density