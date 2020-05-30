from PPV import *
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import pandas as pd


class Root(tk.Tk):
    HEIGHT = 50
    WIDTH = 570

    def __init__(self):
        super(Root, self).__init__()
        canvas = tk.Canvas(self, height=self.HEIGHT)

        tabControl = ttk.Notebook(self)
        tab1 = ttk.Frame(tabControl)
        tabControl.add(tab1, text="Material properties")
        tab2 = ttk.Frame(tabControl)
        tabControl.add(tab2, text="Simulation parameters")

        self.material_panel(tab1)
        self.simulation_panel(tab2)
        tabControl.pack(expan=1, fill="both")

        button = tk.Button(self, text="Run", command=lambda: self.start_simulation())
        button.place(relx=0.35, rely=0.915, relwidth=0.3)
        canvas.pack()

    def material_panel(self, tab):
        self.nshlabel = tk.Label(tab, text="N_Sh: ")
        self.nshentry = tk.Entry(tab)
        self.nshentry.insert(tk.END, "6e17")
        self.ntlabel = tk.Label(tab, text="N_T: ")
        self.ntentry = tk.Entry(tab)
        self.ntentry.insert(tk.END, "5e17")
        self.uphlabel = tk.Label(tab, text="u_ph: ")
        self.uphentry = tk.Entry(tab)
        self.uphentry.insert(tk.END, "5e12")
        self.selabel = tk.Label(tab, text="sigma e: ")
        self.seentry = tk.Entry(tab)
        self.seentry.insert(tk.END, "1e-16")
        self.shlabel = tk.Label(tab, text="sigma h: ")
        self.shentry = tk.Entry(tab)
        self.shentry.insert(tk.END, "1e-15")
        self.ehclabel = tk.Label(tab, text="E_HC: ")
        self.ehcentry = tk.Entry(tab)
        self.ehcentry.insert(tk.END, "0.35")
        self.ehelabel = tk.Label(tab, text="E_HE: ")
        self.eheentry = tk.Entry(tab)
        self.eheentry.insert(tk.END, "0.85")
        self.eeclabel = tk.Label(tab, text="E_EC: ")
        self.eecentry = tk.Entry(tab)
        self.eecentry.insert(tk.END, "0.1")

        self.nshlabel.grid(row=0, column=0, padx=80, pady=15)
        self.nshentry.grid(row=0, column=1, padx=25, pady=15)
        self.ntlabel.grid(row=1, column=0, padx=15, pady=15)
        self.ntentry.grid(row=1, column=1, padx=15, pady=15)
        self.uphlabel.grid(row=2, column=0, padx=15, pady=15)
        self.uphentry.grid(row=2, column=1, padx=15, pady=15)
        self.selabel.grid(row=3, column=0, padx=15, pady=15)
        self.seentry.grid(row=3, column=1, padx=15, pady=15)
        self.shlabel.grid(row=4, column=0, padx=15, pady=15)
        self.shentry.grid(row=4, column=1, padx=15, pady=15)
        self.ehclabel.grid(row=5, column=0, padx=15, pady=15)
        self.ehcentry.grid(row=5, column=1, padx=15, pady=15)
        self.ehelabel.grid(row=6, column=0, padx=15, pady=15)
        self.eheentry.grid(row=6, column=1, padx=15, pady=15)
        self.eeclabel.grid(row=7, column=0, padx=15, pady=15)
        self.eecentry.grid(row=7, column=1, padx=15, pady=15)

    def simulation_panel(self, tab):
        self.nlabel = tk.Label(tab, text="n: ")
        self.nentry = tk.Entry(tab)
        self.nentry.insert(tk.END, "1e16")
        self.dtlabel = tk.Label(tab, text="dt: ")
        self.dtentry = tk.Entry(tab)
        self.dtentry.insert(tk.END, "1")
        self.lstlabel = tk.Label(tab, text="LS time: ")
        self.lstentry = tk.Entry(tab)
        self.lstentry.insert(tk.END, "60")
        self.tminlabel = tk.Label(tab, text="T init: ")
        self.tminentry = tk.Entry(tab)
        self.tminentry.insert(tk.END, "200")
        self.tmaxlabel = tk.Label(tab, text="T final: ")
        self.tmaxentry = tk.Entry(tab)
        self.tmaxentry.insert(tk.END, "250")
        self.mlabel = tk.Label(tab, text="Number of steps: ")
        self.mentry = tk.Entry(tab)
        self.mentry.insert(tk.END, "2")
        self.coollabel = tk.Label(tab, text="Fast cooling")
        self.fastcool = tk.IntVar()
        self.coolbutton = tk.Checkbutton(tab, variable=self.fastcool, command=lambda: self.change_teq_entry_state())
        self.coolbutton.select()
        self.teqlabel = tk.Label(tab, text="Equilibrium temperature: ")
        self.teqentry = tk.Entry(tab)
        self.teqentry.insert(tk.END, "300")
        self.steadystatelabel = tk.Label(tab, text="Compare with steady state")
        self.steadystate = tk.IntVar()
        self.steadystatebutton = tk.Checkbutton(tab, variable=self.steadystate)

        self.nlabel.grid(row=0, column=0, padx=80, pady=15)
        self.nentry.grid(row=0, column=1, padx=25, pady=15)
        self.dtlabel.grid(row=1, column=0, padx=15, pady=15)
        self.dtentry.grid(row=1, column=1, padx=15, pady=15)
        self.lstlabel.grid(row=2, column=0, padx=15, pady=15)
        self.lstentry.grid(row=2, column=1, padx=15, pady=15)
        self.tminlabel.grid(row=3, column=0, padx=15, pady=15)
        self.tminentry.grid(row=3, column=1, padx=15, pady=15)
        self.tmaxlabel.grid(row=4, column=0, padx=15, pady=15)
        self.tmaxentry.grid(row=4, column=1, padx=15, pady=15)
        self.mlabel.grid(row=5, column=0, padx=15, pady=15)
        self.mentry.grid(row=5, column=1, padx=15, pady=15)
        self.coollabel.grid(row=6, column=0, padx=50, pady=15)
        self.coolbutton.grid(row=6, column=1, padx=15, pady=15)
        self.teqlabel.grid(row=7, column=0, padx=50, pady=15)
        self.teqentry.grid(row=7, column=1, padx=15, pady=15)
        self.steadystatelabel.grid(row=8, column=0, padx=50, pady=15)
        self.steadystatebutton.grid(row=8, column=1, padx=15, pady=15)

    def change_teq_entry_state(self):
        if self.teqentry['state'] == tk.NORMAL:
            self.teqentry.configure(state=tk.DISABLED)
        else:
            self.teqentry.configure(state=tk.NORMAL)

    def start_simulation(self):
        ppv = PPV(light=float(self.nentry.get()), temp=300)
        ppv.N_T = float(self.ntentry.get())
        ppv.N_Sh = float(self.nshentry.get())
        ppv.u_ph = float(self.uphentry.get())
        ppv.sigma_h = float(self.shentry.get())
        ppv.sigma_e = float(self.seentry.get())
        ppv.E_HC = float(self.ehcentry.get())
        ppv.E_HE = float(self.eheentry.get())
        ppv.E_EC = float(self.eecentry.get())

        if self.fastcool.get():
            cooltemp = int(self.teqentry.get())
        else:
            cooltemp = None

        results = ppv.ls_wide_simulation(Ti=float(self.tminentry.get()),
                                         Tf=float(self.tmaxentry.get()), m=int(self.mentry.get()),
                                         T0=cooltemp, lstime=float(self.lstentry.get()))

        temperatures = np.linspace(int(self.tminentry.get()), int(self.tmaxentry.get()), int(self.mentry.get()))
        finaloc = [elem[-1] for elem in results]

        plt.subplot(1, 2, 1)

        if self.steadystate.get():
            steadystates = [ppv.eq_state(t) for t in temperatures]
            plt.scatter(temperatures, steadystates, label='steady state')

        plt.scatter(temperatures, finaloc, label='after '+self.lstentry.get()+' s LS')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=2, mode="expand", borderaxespad=0.)
        plt.yscale("log")
        plt.xlabel("T (K)")
        plt.ylabel("p (cm-3)")
        plt.subplot(1, 2, 2)
        plt.xlabel("t (s)")
        plt.ylabel("p (cm-3)")

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'lime', 'orange']

        for table, temp, col in zip(results, temperatures, colors):
            plt.plot(np.linspace(0, float(self.lstentry.get()), np.size(table)), table, col, label=str(int(temp))+'K')

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)

        # axcut = plt.axes([0.7, 0.9, 0.2, 0.075])
        # savebutton = Button(axcut, "Export to CSV")
        # savebutton.on_clicked(lambda x: Root.save_file(np.linspace(0, float(self.lstentry.get()), np.size(results)), results, temperatures))
        plt.show()

    @staticmethod
    def save_file(x, y, temps):
        df = pd.DataFrame(np.transpose(y)[:-1], index=x, columns=temps)
        export_file_path = tk.filedialog.asksaveasfilename(defaultextension='.csv')
        df.to_csv(export_file_path, index=True, header=False)
