import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel, QTabWidget
import Bond_fast
import speciation_fast
import speciation_and_angles
import VaspParser_ML
import LAMMPSParser_Standard
import vibr_spectrum_umd
import umd_to_lammps
import gofr_umd
import msd_umd_fast
import viscosity_new
from PyQt5.QtWidgets import QCheckBox, QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget, QFileDialog, QRadioButton, QButtonGroup, QSpacerItem,QHBoxLayout
import os


def isfloat(n):
    try :
        float(n)
        return True
    except ValueError :
        return False

def usefunction(funct,argv,message):
    try :
        funct.main(argv)
        return True
    except Exception as e :
        message.setText("The error : < "+str(e)+" > has occured. Please check the validity of the file and arguments you provided.")
        return False


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Magmatix : a tool for the use of raw VASP and LAMMPS data files")
        self.setMinimumSize(0,0)
        self.resize(10,10)
        
        
        mainWidget = QWidget()        
        mainLayout = QVBoxLayout()
        main_tab = QTabWidget()
        mainLayout.addWidget(main_tab)
        mainWidget.setLayout(mainLayout)

        # Création du widget à onglets
        
        convert_widget=QWidget()
        convert_layout=QVBoxLayout()
        convert_tab=QTabWidget()
        convert_layout.addWidget(convert_tab)
        convert_widget.setLayout(convert_layout)
        
        
        operations_widget=QWidget()
        operations_layout=QVBoxLayout()
        operations_tab=QTabWidget()
        operations_layout.addWidget(operations_tab)
        operations_widget.setLayout(operations_layout)
        
        calculations_widget=QWidget()
        calculations_layout=QVBoxLayout()
        calculations_tab=QTabWidget()
        calculations_layout.addWidget(calculations_tab)
        calculations_widget.setLayout(calculations_layout)

        
        
        # Création des onglets
        UMD_widget = QWidget()
        Bond_widget = QWidget()
        Speciation_widget = QWidget()
        LAMMPSParser_widget = QWidget()
        Vibration_widget = QWidget()
        viscosity_widget = QWidget()        
        msd_widget = QWidget()
        UMDtoLAMMPS_widget = QWidget()
        gofr_widget = QWidget()
        msd_widget = QWidget()        


        # Ajout des onglets au widget à onglets
        convert_tab.addTab(UMD_widget, "VASP to UMD")
        convert_tab.addTab(LAMMPSParser_widget, "LAMMPS to UMD")
        convert_tab.addTab(UMDtoLAMMPS_widget, "UMD to LAMMPS")
        operations_tab.addTab(gofr_widget,"Creation of the average bonds file")
        operations_tab.addTab(Bond_widget, "Creation of the bonding file")
        operations_tab.addTab(Speciation_widget, "Creation of the population file")
        calculations_tab.addTab(Vibration_widget,"Vibrational spectrum")
        calculations_tab.addTab(viscosity_widget,"Viscosity parameters")
        calculations_tab.addTab(msd_widget,"Mean square deviation")
        

        # Création du contenu des onglets
        self.create_Vibration_layout(Vibration_widget)
        #self.create_Viscosity_layout(Viscosity_widget)
        #self.create_msd_layout(msd_widget)
        self.create_UMD_layout(UMD_widget)
        self.create_LAMMPS_layout(LAMMPSParser_widget)
        self.create_Bond_layout(Bond_widget)
        self.create_Speciation_layout(Speciation_widget)
        self.create_UMD_to_LAMMPS_layout(UMDtoLAMMPS_widget)
        self.create_gofr_layout(gofr_widget)
        self.create_msd_layout(msd_widget)
        self.create_viscosity_layout(viscosity_widget)
        # Définition du widget à onglets comme widget central de la fenêtre principale
        main_tab.addTab(convert_widget,"Conversions")
        main_tab.addTab(operations_widget,"Operations")
        main_tab.addTab(calculations_widget,"Calculations")
        
        self.setCentralWidget(main_tab)
        
    def create_viscosity_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfile_vsc_edit = QLineEdit()
        self.umdfile_vsc=""
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.i_vsc = QLineEdit("0")
        self.n_vsc = QLineEdit("2000")
        
        selectButton.clicked.connect(self.select_UMDfile_vsc)
        computeButton.clicked.connect(self.compute_vsc)
        self.umdfile_vsc_edit.textChanged.connect(self.changeumd_vsc)
        self.vscmessage=QLabel("")
        
        hboxIS = QHBoxLayout()
        hboxIS.addWidget(QLabel("Enter the index of the initial step :"))
        hboxIS.addWidget(self.i_vsc)
        
        hboxLA = QHBoxLayout()
        hboxLA.addWidget(QLabel("Enter the value of lag for autocorrelation (fs ; defalut 2000) :"))
        hboxLA.addWidget(self.n_vsc)
        layout.addWidget(self.umdfile_vsc_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,70))
        layout.addLayout(hboxIS)
        layout.addItem(QSpacerItem(0,70))
        layout.addLayout(hboxLA)
        layout.addItem(QSpacerItem(0,70))
        layout.addWidget(computeButton)
        layout.addWidget(self.vscmessage)
        
        tab_widget.setLayout(layout)
        
        
    def select_UMDfile_vsc(self):
        file_dialog = QFileDialog(self)
        self.umdfile_vsc, _ = file_dialog.getOpenFileName()
        self.umdfile_vsc_edit.setText(self.umdfile_vsc) 

    def changeumd_vsc(self):
        self.umdfile_vsc = self.umdfile_vsc_edit.text()

    def compute_vsc(self):
        if not os.path.isfile(self.umdfile_vsc):
            self.vscmessage.setText("Error : the UMD file "+self.umdfile_vsc+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.i_vsc.text().isnumeric() or int(self.n_vsc.text())<0:
            self.vscmessage.setText("Error : the index of the initial step must be a positive integer.")
        elif not self.n_vsc.text().isnumeric() or int(self.n_vsc.text())<0:
            self.vscmessage.setText("Error : the vertical step must be a positive integer.")
        else :
            argv = ["-f",self.umdfile_vsc,"-i",self.i_vsc.text(),"-n",self.n_vsc.text()]
            result = usefunction(viscosity_new,argv,self.vscmessage)
            if result==True :
                Name = self.umdfile_vsc.split("/")
                name = Name[-1][:-8]+".visc.dat"
                self.vscmessage.setText("Viscosity file successfully created under the name "+name)

        
        
    def create_msd_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfile_msd_edit = QLineEdit()
        self.umdfile_msd=""
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.z_msd = QLineEdit("1")
        self.v_msd = QLineEdit("1")
        self.b_msd = QLineEdit("0")
        
        selectButton.clicked.connect(self.select_UMDfile_msd)
        computeButton.clicked.connect(self.compute_msd)
        self.umdfile_msd_edit.textChanged.connect(self.changeumd_msd)
        self.msdmessage=QLabel("")
        self.atoms = QRadioButton("msd of individual atoms")
        self.elements = QRadioButton("msd of each element")
        buttongroup = QButtonGroup()
        buttongroup.addButton(self.atoms)
        buttongroup.addButton(self.elements)
        self.atoms.setChecked(True)
        
        hboxRB = QHBoxLayout()
        hboxRB.addWidget(self.elements)
        hboxRB.addWidget(self.atoms)
        
        hboxHJ = QHBoxLayout()
        hboxHJ.addWidget(QLabel("Enter the value of the horizontal jump :"))
        hboxHJ.addWidget(self.z_msd)

        hboxVJ = QHBoxLayout()
        hboxVJ.addWidget(QLabel("Enter the value of the vertical jump :"))
        hboxVJ.addWidget(self.v_msd)
        
        hboxBR = QHBoxLayout()
        hboxBR.addWidget(QLabel("Enter the lenght of the ballistic regime :"))
        hboxBR.addWidget(self.b_msd)

        layout.addWidget(self.umdfile_msd_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxHJ)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxVJ)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxBR)
        layout.addItem(QSpacerItem(0,-10))
        layout.addWidget(QLabel("Choose the form of the data to be displayed in the output file :"))
        layout.addItem(QSpacerItem(0,-20))
        layout.addLayout(hboxRB)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(computeButton)
        layout.addItem(QSpacerItem(0,-20))
        layout.addWidget(self.msdmessage)
        
        tab_widget.setLayout(layout)
    
        
    def select_UMDfile_msd(self):
        file_dialog = QFileDialog(self)
        self.umdfile_msd, _ = file_dialog.getOpenFileName()
        self.umdfile_msd_edit.setText(self.umdfile_msd) 

    def changeumd_msd(self):
        self.umdfile_msd = self.umdfile_msd_edit.text()

    def compute_msd(self):
        if not os.path.isfile(self.umdfile_msd):
            self.msdmessage.setText("Error : the UMD file "+self.umdfile_msd+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.z_msd.text().isnumeric() or int(self.z_msd.text())<1:
            self.msdmessage.setText("Error : the horizontal jump must be a strictly positive integer.")
        elif not self.v_msd.text().isnumeric() or int(self.v_msd.text())<1:
            self.msdmessage.setText("Error : the vertical step must be a strictly positive integer.")
        elif not self.b_msd.text().isnumeric() or int(self.b_msd.text())<0:
            self.msdmessage.setText("Error : the length of the ballistic regime must be a positive integer.")
        else :
            if self.atoms.isChecked():                
                argv = ["-f",self.umdfile_msd,"-z",self.z_msd.text(),"-v",self.v_msd.text(),"-b",self.b_msd.text(),"-x","atoms"]
            else :
                argv = ["-f",self.umdfile_msd,"-z",self.z_msd.text(),"-v",self.v_msd.text(),"-b",self.b_msd.text(),"-x","elements"]
            result = usefunction(msd_umd_fast,argv,self.msdmessage)
            if result==True :
                Name = self.umdfile_msd.split("/")
                name = Name[-1][:-8]+".msd.dat"
                self.msdmessage.setText("msd file successfully created under the name "+name)
        
    
    
    
    
    
    
    
        

    def create_gofr_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfileg_edit = QLineEdit()
        self.umdfileg=""
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.s_g = QLineEdit("1")
        self.d_g = QLineEdit("1")
        self.i_g = QLineEdit("0")
        self.gpu = False

        self.gpubox = QCheckBox("Use gpu portage (be sure your computer is configured to handle it !)")
        self.gpubox.stateChanged.connect(self.changedgpubox)
        
        selectButton.clicked.connect(self.select_UMDfile_gofr)
        computeButton.clicked.connect(self.compute_gofr)
        self.umdfileg_edit.textChanged.connect(self.changeumdg)
        self.gofrmessage=QLabel("")
        
        hboxSF = QHBoxLayout()
        hboxSF.addWidget(QLabel("Enter the sampling frequency :"))
        hboxSF.addWidget(self.s_g)
        
        hboxDI = QHBoxLayout()
        hboxDI.addWidget(QLabel("Enter the discretization interval :"))
        hboxDI.addWidget(self.d_g)
        
        hboxIS = QHBoxLayout()
        hboxIS.addWidget(QLabel("Enter the index of the initial step :"))
        hboxIS.addWidget(self.i_g)
        
        layout.addWidget(self.umdfileg_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxSF)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxDI)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxIS)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.gpubox)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(computeButton)
#        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.gofrmessage)
        
        tab_widget.setLayout(layout)
        
        
    def select_UMDfile_gofr(self):
        file_dialog = QFileDialog(self)
        self.umdfileg, _ = file_dialog.getOpenFileName()
        self.umdfileg_edit.setText(self.umdfileg) 

    def changeumdg(self):
        self.umdfileg = self.umdfileg_edit.text()

    def changedgpubox(self,state):
        if state == 2:
            self.gpu = True
        else : 
            self.gpu = False

    def compute_gofr(self):
        if not os.path.isfile(self.umdfileg):
            self.gofrmessage.setText("Error : the file "+self.umdfileg+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.s_g.text().isnumeric() or int(self.s_g.text())<1:
            self.gofrmessage.setText("Error : the sampling frequency must be a strictly positive integer.")
        elif not self.d_g.text().isnumeric() or int(self.d_g.text())<1:
            self.gofrmessage.setText("Error : the discretization interval must be a strictly positive integer.")
        elif not self.i_g.text().isnumeric() or int(self.i_g.text())<0:
            self.gofrmessage.setText("Error : the initial step must be a positive integer.")
        else :
            argv = ["-f",self.umdfileg,"-s",self.s_g.text(),"-d",self.d_g.text(),"-i",self.i_g.text(),"-g","use_gpu="+str(self.gpu)]
            result = usefunction(gofr_umd,argv,self.gofrmessage)
            if result==True :
                Name = self.umdfileg.split("/")
                name = Name[-1][:-8]+".gofr.dat"
                self.gofrmessage.setText("gofr file successfully created under the name "+name)








    def create_Vibration_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.UMDfileVib=""
        self.messageVib=QLabel("")
        self.temperature = "5000"
        self.temperature_edit = QLineEdit("5000")
        self.temperature_edit.textChanged.connect(self.changetemp)
        self.UMDfileVib_edit = QLineEdit()
        
        UMDButton = QPushButton("Select the UMD file")
        
        UMDButton.clicked.connect(self.select_UMDfileVib)
        self.UMDfileVib_edit.textChanged.connect(self.changevib)

        hboxT=QHBoxLayout()
        hboxT.addWidget(QLabel("Select the temperature (default 5000 K)"))
        hboxT.addWidget(self.temperature_edit)
                
        
        ComputeVibButton = QPushButton("Compute vibrational spectrum")
        ComputeVibButton.clicked.connect(self.ComputeVib)

        layout.addWidget(self.UMDfileVib_edit)
        layout.addWidget(UMDButton)
        layout.addItem(QSpacerItem(0,70))
        layout.addLayout(hboxT)
        layout.addItem(QSpacerItem(0,70))
        layout.addWidget(ComputeVibButton)
        layout.addItem(QSpacerItem(0,-10))
        layout.addWidget(self.messageVib)
        tab_widget.setLayout(layout)

    def ComputeVib(self):
        argv = ["-f",self.UMDfileVib,"-t",self.temperature]
        
        if self.UMDfileVib=="":
            self.messageVib.setText("Error : Please select an UMD file.")
        elif not os.path.isfile(self.UMDfileVib_edit.text()):
            self.messageVib.setText("Error : the file "+self.UMDfileVib+" is displaced or missing.") 
        elif not (self.temperature.isnumeric() and int(self.temperature)>=0):
            self.messageVib.setText("Error : the temperature must be a strictly positive integer.")
        else :
            self.messageVib.setText("Computing...") 
            result=usefunction(vibr_spectrum_umd.main,argv,self.messageVib)
            if result==True :
                self.messageVib.setText("Vibrational spectrum successfully calculated. Files created under the names "+self.UMDfileVib[:-8]+".vels.scf.dat and "+self.UMDfileVib[:-8]+".vibr.dat")

    def select_UMDfileVib(self):
        file_dialog = QFileDialog(self)
        self.UMDfileVib, _ = file_dialog.getOpenFileName()
        self.UMDfileVib_edit.setText(self.UMDfileVib) 
    
    def changetemp(self):
        self.temperature = self.temperature_edit.text()

    def changevib(self):
        self.UMDfileVib = self.UMDfileVib_edit.text()









    def create_LAMMPS_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.LAMMPSfile=""
        self.logfile=""
        self.pressfile=""
        
        self.LAMMPSfile_edit = QLineEdit()
        self.logfile_edit = QLineEdit()
        self.pressfile_edit = QLineEdit()
        
        
        self.LAMMPSfile_edit.textChanged.connect(self.changeLAMMPS)
        self.logfile_edit.textChanged.connect(self.changelog)
        self.pressfile_edit.textChanged.connect(self.changepress)
        
        
        selectLAMMPSButton = QPushButton("Select the LAMMPS file")
        selectLAMMPSButton.clicked.connect(self.select_LAMMPS)
        selectlogButton = QPushButton("Select the log file")
        selectlogButton.clicked.connect(self.select_log)
        selectpressButton = QPushButton("Select the press file")
        selectpressButton.clicked.connect(self.select_press)
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_LAMMPS)
        
        self.L2Umessage = QLabel("")
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.LAMMPSfile_edit)
        layout.addWidget(selectLAMMPSButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.logfile_edit)
        layout.addWidget(selectlogButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.pressfile_edit)
        layout.addWidget(selectpressButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.L2Umessage)
        
        tab_widget.setLayout(layout)


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
        if not os.path.isfile(self.LAMMPSfile):
            self.L2Umessage.setText("ERROR : LAMMPS file not found. Please check the path or the file name.")        
        elif not os.path.isfile(self.logfile):
            self.L2Umessage.setText("ERROR : log file not found. Please check the path or the file name.")        
        elif not os.path.isfile(self.pressfile):
            self.L2Umessage.setText("ERROR : press file not found. Please check the path or the file name.")        
        else:            
            argv=["-f",self.LAMMPSfile,"-l",self.logfile,"-a",self.pressfile]
            result = usefunction(LAMMPSParser_Standard,argv,self.L2Umessage)
            if result :
                Name = self.LAMMPSfile.split("/")
                name = Name[-1]+".umd.dat"
                self.L2Umessage.setText("The UMD file has been successfully created under the name "+name)

    def changeLAMMPS(self):
        self.LAMMPSfile = self.LAMMPSfile_edit.text()
    def changelog(self):
        self.logfile = self.logfile_edit.text()
    def changepress(self):
        self.pressfile = self.pressfile_edit.text()







    def create_UMD_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.VASPfile_edit = QLineEdit()
        self.VASPfile = ""
        
        selectVASPButton = QPushButton("Select the VASP file")
        selectVASPButton.clicked.connect(self.select_VASP)
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_VASP)
        
        self.VASPmessage = QLabel("")
        
        self.VASPfile_edit.textChanged.connect(self.changeVASP)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.VASPfile_edit)
        layout.addWidget(selectVASPButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.VASPmessage)
        
        tab_widget.setLayout(layout)
        
    def select_VASP(self):
        file_dialog = QFileDialog(self)
        self.VASPfile, _ =file_dialog.getOpenFileName()
        self.VASPfile_edit.setText(self.VASPfile)
        
    
    def convert_VASP(self):
        if not os.path.isfile(self.VASPfile):
            self.VASPmessage.setText("ERROR : the VASP file "+self.VASPfile+ " is displaced or missing. Please check the path or filename.")        
        else:
            argv=["-f",self.VASPfile]
            result = usefunction(VaspParser_ML,argv,self.VASpmessage)
            if result :
                Name = self.VASPfile.split("/")
                name = Name[-1]+".umd.dat"
                self.VASPmessage.setText("The UMD file has been successfully created under the name "+name)


    def changeVASP(self):
        self.VASPfile = self.VASPfile_edit.text()   




            
    def create_UMD_to_LAMMPS_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.U2Lfile_edit = QLineEdit()
        self.U2L = ""
        
        selectUMDButton = QPushButton("Select the UMD file")
        selectUMDButton.clicked.connect(self.select_U2L)
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_U2L)
        
        self.U2Lmessage = QLabel("")
        
        self.U2Lfile_edit.textChanged.connect(self.changeU2L)

        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.U2Lfile_edit)
        layout.addWidget(selectUMDButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.U2Lmessage)
        
        tab_widget.setLayout(layout)
        
    def select_U2L(self):
        file_dialog = QFileDialog(self)
        self.U2Lfile, _ =file_dialog.getOpenFileName()
        self.U2Lfile_edit.setText(self.U2Lfile)
        
    
    def convert_U2L(self):
        if not os.path.isfile(self.U2Lfile):
            self.U2Lmessage.setText("ERROR : UMD file "+self.U2Lfile+" is displaced or missing. Please check the path or the file name.")        
        else:
            argv=["-f",self.U2Lfile]
            result = usefunction(umd_to_lammps,argv,self.U2LMessage)
            if result :
                Name = self.U2Lfile.split("/")
                name = Name[-1][:-8]+".lammps"
                self.U2Lmessage.setText("The LAMMPS file has been successfully created under the name "+name)

    def changeU2L(self):
        self.U2Lfile = self.U2Lfile_edit.text()




    

    def create_Speciation_layout(self, tab_widget):
        layout = QVBoxLayout()
        
        self.bondfile_edit = QLineEdit()
        self.bondfile =""
        self.centralatom_edit = QLineEdit()
        self.outeratom_edit = QLineEdit()
        self.umdfileSpec=""
        self.umdfileSpec_edit = QLineEdit()
        self.umdfileSpec_edit.hide()
        self.rings = QCheckBox("All levels of coordination")
        self.angles = QCheckBox("Calculate the angles whithin each polyhedra (needs to have, and sets, a degree of coordination equal to 1)")
        self.rings_edit = QLineEdit()
        
        self.ring = 1
        self.angle = False
        
        
        self.bondfile_edit.textChanged.connect(self.changebonds)

        selectbondbutton = QPushButton("Select the bond file")
        selectbondbutton.clicked.connect(self.select_bondfile)
        self.specmessage = QLabel("")
        
        submit_button = QPushButton("Compute")
        submit_button.clicked.connect(self.compute_speciation)
        
        self.rings.stateChanged.connect(self.changering)
        self.angles.stateChanged.connect(self.changeangles)
        
        self.rings_edit.textChanged.connect(self.check_boxes)
        
        self.umdfileSpec_edit.textChanged.connect(self.changeumdSpec)
        self.umdSpecbutton = QPushButton("Select the matching umd file to compute the angles.")
        self.umdSpecbutton.hide()
        self.umdSpecbutton.clicked.connect(self.select_umdfileSpec)
        
        hboxCA=QHBoxLayout()
        hboxCA.addWidget(QLabel("Central atoms element :"))
        hboxCA.addWidget(self.centralatom_edit)
        
        hboxOA=QHBoxLayout()
        hboxOA.addWidget(QLabel("Outer atoms element :"))
        hboxOA.addWidget(self.outeratom_edit)

        
        layout.addWidget(self.bondfile_edit)
        layout.addWidget(selectbondbutton)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxCA)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxOA)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(QLabel("Enter the degree of coordination :"))
        layout.addWidget(self.rings_edit)
        layout.addWidget(self.rings)
        layout.addWidget(self.angles)
        layout.addItem(QSpacerItem(0,15))
        layout.addWidget(self.umdfileSpec_edit)
        layout.addWidget(self.umdSpecbutton)
        layout.addItem(QSpacerItem(0,30))
        layout.addWidget(submit_button)
        layout.addItem(QSpacerItem(0,30))
        layout.addWidget(self.specmessage)
        
#        layout.setSpacing(2)
        
        tab_widget.setLayout(layout)

    def check_boxes(self):
        if self.rings_edit.text()!="1" and self.angle:
            self.angles.toggle()
        if self.rings_edit.text()!="0" and self.rings_edit.text()!="All levels (0)" and self.rings.isChecked() :
            self.rings.toggle()

    def changering(self,state):
        if state == 2:
            self.ring = 0
            self.rings_edit.setText("All levels (0)")
        else :
            if not self.rings_edit.text().isnumeric():
                self.ring = 1
                self.rings_edit.setText("1")
                
    def changeangles(self,state):
        if state == 2:
            self.ring = 1
            self.angle = True
            self.rings_edit.setText("1")
            self.umdfileSpec_edit.show()
            self.umdSpecbutton.show()
        else :
            self.angle = False
            self.umdfileSpec_edit.hide()
            self.umdSpecbutton.hide()

    def select_bondfile(self):
        file_dialog = QFileDialog(self)
        self.bondfile, _ = file_dialog.getOpenFileName()        
        self.bondfile_edit.setText(self.bondfile)

    def select_umdfileSpec(self):
        file_dialog = QFileDialog(self)
        self.umdfileSpec, _ = file_dialog.getOpenFileName()        
        self.umdfileSpec_edit.setText(self.umdfileSpec)    

    def compute_speciation(self):
        
        if self.rings_edit.text().isnumeric():
            self.ring = self.rings_edit.text()
        elif self.rings_edit.text() != "All levels (0)":
            self.rings_edit.setText("1")
            self.ring=1
        
        if os.path.isfile(self.bondfile_edit.text()):
            self.specmessage.setText("")
            if self.angle :
                if os.path.isfile(self.umdfileSpec):                    
                    argv=["-f",self.bondfile,"-u",self.umdfileSpec,"-c",self.centralatom_edit.text(),"-a",self.outeratom_edit.text(),"-r", self.ring]
                    speciation_and_angles.main(argv)
                elif self.umdfileSpec=="" :
                    self.specmessage.setText("If you want to calculate the angles, please select an UMD file before launching computation.")            
                else :
                    self.specmessage.setText("ERROR : the file "+self.umdfileSPec+" is displaced or missing")
            else :
                if self.outeratom_edit.text()=="":
                    self.specmessage.setText("Please select an element for the outer atom.")
                elif self.centralatom_edit.text()=="":
                    self.specmessage.setText("Please select an element for the central atom.")
                else :
                    argv=["-f",self.bondfile,"-c",self.centralatom_edit.text(),"-a",self.outeratom_edit.text(),"-r", self.ring]
                    result=usefunction(speciation_fast,argv,self.specmessage)
                    if result :
                        self.specmessage.setText("Population file successfully created under the name "+self.bondfile_edit.text().split("/")[-1][:-4]+".r"+str(self.ring)+".popul.dat")
        elif self.bondfile_edit.text()=="":
            self.specmessage.setText("Please select a bond file before launching computation.")            
        else:
            self.specmessage.setText("ERROR : the file "+self.bondfile_edit.text()+" is displaced or missing")

    def changebonds(self):
        self.bondfile = self.bondfile_edit.text()

    def changeumdSpec(self):
        self.umdfileSpec = self.umdfileSpec_edit.text()








    def create_Bond_layout(self, tab_widget):
        layout = QVBoxLayout()
        # Création des widgets de saisie
        self.umdfile_edit = QLineEdit()
        self.inpfile_edit = QLineEdit()
        
        self.umdfile=""
        self.inpfile=""
        
        self.umdfile_edit.textChanged.connect(self.changeumd)
        self.inpfile_edit.textChanged.connect(self.changeinp)
                
        self.bondmessage = QLabel("")
        
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

        hboxBL = QHBoxLayout()
        hboxBL.addWidget(QLabel("Bonding length in Angströms (optional ; will overwrite the values of the input file) :"))
        hboxBL.addWidget(self.l_edit)

        hboxSF = QHBoxLayout()        
        hboxSF.addWidget(QLabel("Sampling frequency :"))
        hboxSF.addWidget(self.s_edit)

        # Création d'un layout vertical pour organiser les widgets
        layout.addWidget(self.umdfile_edit)
        layout.addWidget(select_button)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.inpfile_edit)
        layout.addWidget(select_buttonI)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxBL)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxSF)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.bondmessage)

        # Création d'un widget de base et définition du layout comme layout principal
        widget = QWidget()
        widget.setLayout(layout)

        # Définition du widget comme widget central de la fenêtre principale
        self.setCentralWidget(widget)
        
        tab_widget.setLayout(layout)
        
    def compute_Bond(self):
           
        if not os.path.isfile(self.umdfile):
            if self.umdfile=="":
                self.bondmessage.setText("Please select an umd file before launching computation.")
            else :
                self.bondmessage.setText("Error : the file "+self.umdfile_edit.text()+" is displaced or missing. Please check the path or file name.") 
        elif not os.path.isfile(self.inpfile) and not(isfloat(self.l_edit.text())):
            if self.inpfile=="":
                self.bondmessage.setText("Please select an input file or set manually the value of the bonding length before launching the computation.")
            else :
                self.bondmessage.setText("Error : the file "+self.inpfile_edit.text()+" is displaced or missing. Please check the path or file name.") 
        elif not(self.s_edit.text().isnumeric()) or int(self.s_edit.text())<1 :        
            self.bondmessage.setText("The value of the sampling frequency has to be a strictly positive integer.")
        else :
            if not(isfloat(self.l_edit.text())):
                argv=["-f",self.umdfile,"-s",self.s_edit.text(),"-i",self.inpfile]
                self.l_edit=""
                self.bondmessage.setText("The .inp file has been used to determine the bonds.")
            else :
                argv=["-f",self.umdfile,"-s",self.s_edit.text(),"-l",self.l_edit.text(),"-i",self.inpfile]
                self.bondmessage.setText("The explicit value of "+self.l_edit.text()+" (Angströms) has been used to determine the bonds.")

            self.bondmessage.setText("Computing...")
            result = usefunction(Bond_fast,argv,self.bondname)
            if result :
                self.bondmessage.setText("Bonds file successfully created under the name "+self.umdfile.split("/")[-1][:-8]+".bondingfile.dat") 
                self.bondfile_edit.setText(self.umdfile[:-8]+".bondingfile.dat")


    def select_inpfile(self):
        file_dialog = QFileDialog(self)
        self.inpfile, _ = file_dialog.getOpenFileName()
        self.inpfile_edit.setText(self.inpfile)
        
    def select_UMDfile(self):
        file_dialog = QFileDialog(self)
        self.umdfile, _ = file_dialog.getOpenFileName()
        self.umdfile_edit.setText(self.umdfile)
    
    def changeinp(self):
        self.inpfile = self.inpfile_edit.text()
    def changeumd(self):
        self.umdfile = self.umdfile_edit.text()
    
    
    
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
