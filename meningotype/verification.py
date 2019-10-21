import pathlib, pandas, subprocess, datetime,os, numpy


class Verification:

    def __init__(self,verification_path, reason):

        if pathlib.Path(verification_path).is_absolute():
            self.verification_path = pathlib.Path(verification_path)
        else:
            self.verification_path = pathlib.Path(verification_path).absolute()

        if self.check_verification_path():
            self.assemblies = self.verification_path / 'assemblies'
            self.result_dir = self.verification_path / 'verification_results'
        self.reason = reason
        self.day = datetime.datetime.today().strftime("%d_%m_%y")
        self.reverification_path = self.result_dir / f"reverification_{self.reason}_{self.day}.txt"
        self.original_results = self.verification_path / "reverification_isolates_dummy.txt"

    def check_verification_path(self):
        '''
        check that verification path exists
        '''
        if self.verification_path.exists():
            return True

        else:
            print('Verification path does not exist. Please try again')
            raise SystemExit
    
    def get_sed_cmd(self):
        first_sed = f"{self.assemblies}".replace("/", "\\/")
        if not first_sed.endswith("/"):
            first_sed = f"{first_sed}\\/"
        sed_cmd = f"| sed 's/{first_sed}//g' | sed 's/\.fa//g'"

        return sed_cmd

     
    def run_meningotyper(self):
        '''
        run meningotyper
        '''
        sed_cmd = self.get_sed_cmd()
        cmd = f"python3 meningotype/meningotype.py --all {self.assemblies}/*.fa {sed_cmd} > {self.reverification_path}"
        print(cmd)
        p = subprocess.run(cmd, shell = True, capture_output = True)
        if p.returncode == 0:
            return True
        else:
            print(f'Something the seems to have gone wrong. The following error was returned by meningotype {p.returncode}')
            raise SystemExit
      
    def get_original(self):

        if self.original_results.exists():
            tab = pandas.read_csv(self.original_results, sep = '\t', index_col=False, usecols=['Dummy ID', 'expected_serotype', 'expected_PorA', 'expected_FetA'])
            tab = tab.rename(columns= {"Dummy ID":'SAMPLE_ID'})
            tab['expected_PorA'] = tab['expected_PorA'].apply(lambda x:x.replace('P1.', ''))
            return tab
        else:
            print(f"The original results seem to be missing. File {self.original_results} can not be found. Please make sure you have a file for comparison.")
            raise SystemExit

    def verify(self):
        
        passed =  self.run_meningotyper()

        if passed:
            
            orig = self.get_original()
            # print(self.reverification_path)
            print("Opening meningotype results for file for comparisons")
            reverf = pandas.read_csv(self.reverification_path, sep = '\t', index_col= False, usecols=['SAMPLE_ID', 'SEROGROUP', 'PorA', 'FetA'])
            comp = pandas.merge(orig, reverf)
            print("making comparisons")
            comp['SeroComp'] = numpy.where(comp['SEROGROUP'] == comp['expected_serotype'], 'PASS', 'FAIL')    
            comp['PorAComp'] = numpy.where(comp['PorA'] == comp['expected_PorA'], 'PASS', 'FAIL')    
            comp['FetAComp'] = numpy.where(comp['FetA'] == comp['expected_FetA'], 'PASS', 'FAIL')    
            print("saving comparison files")
            comp.to_csv(self.reverification_path, sep = '\t', index = False)
 