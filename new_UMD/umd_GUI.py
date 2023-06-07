import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel, QTabWidget
import Bond_par_C
import speciation_fast
import VaspParser_ML
import LAMMPSParser_Standard
from PyQt5.QtWidgets import QCheckBox, QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget, QFileDialog
import os

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Magmatix : a tool for the treatment of raw VASP and LAMMPS data files")
        self.resize(600,600)
        
        

        # Création du widget à onglets
        tab_widget = QTabWidget()

        # Création des onglets
        UMD_tab = QWidget()
        Bond_tab = QWidget()
        Speciation_tab = QWidget()
        LAMMPSParser_tab = QWidget()

        # Ajout des onglets au widget à onglets
        tab_widget.addTab(UMD_tab, "VASP to UMD")
        tab_widget.addTab(LAMMPSParser_tab, "LAMMPS to UMD")
        tab_widget.addTab(Bond_tab, "Computation of the bonding file")
        tab_widget.addTab(Speciation_tab, "Computation of the population file")

        # Création du contenu des onglets
        self.create_UMD_tab(UMD_tab)
        self.create_LAMMPS_tab(LAMMPSParser_tab)
        self.create_Bond_tab(Bond_tab)
        self.create_Speciation_tab(Speciation_tab)

        # Définition du widget à onglets comme widget central de la fenêtre principale
        self.setCentralWidget(tab_widget)
        
    def create_LAMMPS_tab(self,tab_widget):
        layout = QVBoxLayout()
        
        self.LAMMPSfile_edit = QLineEdit()
        self.logfile_edit = QLineEdit()
        self.pressfile_edit = QLineEdit()
        
        selectLAMMPSButton = QPushButton("Select the LAMMPS file")
        selectLAMMPSButton.clicked.connect(self.select_LAMMPS)
        selectlogButton = QPushButton("Select the log file")
        selectlogButton.clicked.connect(self.select_log)
        selectpressButton = QPushButton("Select the press file")
        selectpressButton.clicked.connect(self.select_press)
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_LAMMPS)
        
        self.label = QLabel("")
        
        layout.addWidget(self.LAMMPSfile_edit)
        layout.addWidget(selectLAMMPSButton)
        layout.addWidget(self.logfile_edit)
        layout.addWidget(selectlogButton)
        layout.addWidget(self.pressfile_edit)
        layout.addWidget(selectpressButton)
        layout.addWidget(submit_button)
        layout.addWidget(self.label)
        
        tab_widget.setLayout(layout)

    def create_UMD_tab(self,tab_widget):
        layout = QVBoxLayout()
        
        self.VASPfile_edit = QLineEdit()
        
        selectVASPButton = QPushButton("Select the VASP file")
        selectVASPButton.clicked.connect(self.select_VASP)
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_VASP)
        
        self.label = QLabel("")
        
        layout.addWidget(self.VASPfile_edit)
        layout.addWidget(selectVASPButton)
        layout.addWidget(submit_button)
        layout.addWidget(self.label)
        
        tab_widget.setLayout(layout)

    def create_Speciation_tab(self, tab_widget):
        layout = QVBoxLayout()
        
        self.bondfile_edit = QLineEdit()
        self.centralatom_edit = QLineEdit()
        self.outeratom_edit = QLineEdit()
        self.rings = QCheckBox("All levels of coordination")
        self.rings_edit = QLineEdit()
        
        self.ring = 1
        
        
        self.populname = QLabel("")

        selectbondbutton = QPushButton("Select the bond file")
        selectbondbutton.clicked.connect(self.select_bondfile)
        
        
        submit_button = QPushButton("Compute")
        submit_button.clicked.connect(self.compute_speciation)
        
        self.rings.stateChanged.connect(self.changering)
        
        
        layout.addWidget(self.bondfile_edit)
        layout.addWidget(selectbondbutton)
        layout.addWidget(QLabel(""))
        layout.addWidget(QLabel(""))
        layout.addWidget(QLabel("Central atom :"))
        layout.addWidget(self.centralatom_edit)
        layout.addWidget(QLabel("Outer atom :"))
        layout.addWidget(self.outeratom_edit)
        layout.addWidget(QLabel(""))
        layout.addWidget(QLabel("Enter the degree of coordination :"))
        layout.addWidget(self.rings_edit)
        layout.addWidget(self.rings)
        layout.addWidget(QLabel(""))
        layout.addWidget(QLabel(""))
        
        layout.addWidget(submit_button)
        
        layout.addWidget(self.populname)
        layout.setSpacing(2)
        layout.setContentsMargins(0,20,0,200)
        
        tab_widget.setLayout(layout)

    def create_Bond_tab(self, tab_widget):
        layout = QVBoxLayout()
        # Création des widgets de saisie
        self.umdfile_edit = QLineEdit()
        self.inpfile_edit = QLineEdit()
        
        self.bondname = QLabel("")
        
        self.error_message = QLabel("")
        
        self.l_edit = QLineEdit()
        self.s_edit = QLineEdit()
        # Création d'un bouton pour sélectionner un fichier
        select_button = QPushButton("Select the UMD file")
        select_button.clicked.connect(self.select_UMDfile)

        select_buttonI = QPushButton("Select the average bonding input file")
        select_buttonI.clicked.connect(self.select_inpfile)

        
        # Création d'un bouton pour soumettre les paramètres
        submit_button = QPushButton("Compute")
        submit_button.clicked.connect(self.compute_Bond)

        # Création d'un layout vertical pour organiser les widgets
        layout.addWidget(self.umdfile_edit)
        layout.addWidget(select_button)
        layout.addWidget(self.inpfile_edit)
        layout.addWidget(select_buttonI)
        layout.addWidget(QLabel("Bonding length (optional ; will overwrite the values of the input file) :"))
        layout.addWidget(self.l_edit)
        layout.addWidget(QLabel("Sampling frequency :"))
        layout.addWidget(self.s_edit)
        layout.addWidget(submit_button)
        layout.addWidget(self.bondname)
        layout.addWidget(self.error_message)

        # Création d'un widget de base et définition du layout comme layout principal
        widget = QWidget()
        widget.setLayout(layout)

        # Définition du widget comme widget central de la fenêtre principale
        self.setCentralWidget(widget)
        
        tab_widget.setLayout(layout)

    def select_UMDfile(self):
        file_dialog = QFileDialog(self)
        self.umdfile, _ = file_dialog.getOpenFileName()
        self.umdfile_edit.setText(self.umdfile)
        
    def select_inpfile(self):
        file_dialog = QFileDialog(self)
        self.inpfile, _ = file_dialog.getOpenFileName()
        self.inpfile_edit.setText(self.inpfile)

    def compute_Bond(self):
        # Utilisation de la valeur du fichier sélectionné
        if self.umdfile_edit.text()!="" and (self.inpfile_edit.text()!="" or self.l_edit.text()!=""):
            print("Files selected:", self.umdfile_edit.text()," and ",self.inpfile_edit.text())
            argv=["-f",self.umdfile_edit.text(),"-s",self.s_edit.text(),"-l",self.l_edit.text(),"-i",self.inpfile_edit.text()]
            if not os.path.isfile(self.umdfile_edit.text()):
               self.error_message.setText("Error : the file "+self.umdfile_edit.text()+" is displaced or missing.") 
            elif not os.path.isfile(self.inpfile_edit.text()) and self.l_edit.text()=="":
                  self.error_message.setText("Error : the file "+self.inpfile_edit.text()+" is displaced or missing.") 
            else:
                if self.s_edit.text().isnumeric():
                    if int(self.s_edit.text())<1 :
                        self.error_message.setText("The value of the sampling frequency has to be a strictly positive integer.")
                    else :
                        Bond_par_C.main(argv)
                        self.error_message.setText("")
                        self.bondname.setText("Bond file successfully created under the name "+self.umdfile.split("/")[-1][:-8]+".bondingfile.dat") 
                        self.bondfile_edit.setText(self.umdfile[:-8]+".bondingfile.dat")
                else :
                    self.error_message.setText("Please set the value of the sampling frequency with a strictly positive integer.")
        elif self.umdfile_edit.text()!="":
            self.error_message.setText("Please select an input file or set manually the value of the bonding length before launching the computation.")
        else :
            self.error_message.setText("Please select an umd file before launching computation.")


    def select_bondfile(self):
        file_dialog = QFileDialog(self)
        self.bondfile, _ = file_dialog.getOpenFileName()        
        self.bondfile_edit.setText(self.bondfile)

    def changering(self,state):
        if state == 2:
            self.ring = 0
            self.rings_edit.setText("All levels (0)")
        else :
            self.ring = 1
            self.rings_edit.setText("1")
        
    def compute_speciation(self):
        
        if self.rings_edit.text().isnumeric():
            self.ring = self.rings_edit.text()
        elif self.rings_edit.text() != "All levels (0)":
            self.rings_edit.setText("1")
            self.ring=1
        
        if os.path.isfile(self.bondfile_edit.text()):
            argv=["-f",self.bondfile_edit.text(),"-c",self.centralatom_edit.text(),"-a",self.outeratom_edit.text(),"-r", self.ring]
            speciation_fast.main(argv)
            self.populname.setText("Population file successfully created under the name "+self.bondfile_edit.text().split("/")[-1][:-4]+".r"+str(self.ring)+".popul.dat")
        else :
            self.error_message_spec.setText("ERROR : Bonding file not found")
        

        
    def select_LAMMPS(self):
        file_dialog = QFileDialog(self)
        self.LAMMPSfile, _ =file_dialog.getOpenFileName()
        self.LAMMPSfile_edit.setText(self.LAMMPSfile)
        
    def select_log(self):
        file_dialog = QFileDialog(self)
        self.logfile, _ =file_dialog.getOpenFileName()
        self.logfile_edit.setText(self.logfile)

    def select_press(self):
        file_dialog = QFileDialog(self)
        self.pressfile, _ =file_dialog.getOpenFileName()
        self.pressfile_edit.setText(self.pressfile)        
    
    def convert_LAMMPS(self):
        argv=["-f",self.LAMMPSfile,"-l",self.logfile,"-a",self.pressfile]
        LAMMPSParser_Standard.main(argv)
        Name = self.LAMMPSfile.split("/")
        name = Name[-1]+".umd.dat"
        self.label.setText("The UMD file has been successfully created under the name "+name)
        
    def select_VASP(self):
        file_dialog = QFileDialog(self)
        self.VASPfile, _ =file_dialog.getOpenFileName()
        self.VASPfile_edit.setText(self.VASPfile)
        
    
    def convert_VASP(self):
        argv=["-f",self.VASPfile]
        VaspParser_ML.main(argv)
        Name = self.VASPfile.split("/")
        name = Name[-1]+".umd.dat"
        self.label.setText("The UMD file has been successfully created under the name "+name)

        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
