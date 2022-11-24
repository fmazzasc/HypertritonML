import os


# RUN ANALYSIS W/ MATTER + ANTIMATTER 
os.system("python3 run_analysis.py ../Config/2body_analysis_upd_large_bins.yaml -t -a -s")
os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -dbshape -s')
os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -s')
os.system('python3 compute_blambda.py ../Config/2body_analysis_upd_large_bins.yaml -dbshape -s -syst')
os.system('python3 compute_blambda.py ../Config/2body_analysis_upd_large_bins.yaml -s -syst')
os.system('python3 compute_lifetime.py ../Config/2body_analysis_upd_large_bins.yaml -syst -s')

# RUN ANALYSIS W/ SPLITTED MATTER AND ANTIMATTER  
# os.system("python3 run_analysis.py ../Config/2body_analysis_upd_large_bins.yaml -t -a -s --antimatter")
# os.system("python3 run_analysis.py ../Config/2body_analysis_upd_large_bins.yaml -t -a -s --matter")
# os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -s -dbshape --matter')
# os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -s -dbshape --antimatter')
# os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -s --matter')
# os.system('python3 signal_extraction.py ../Config/2body_analysis_upd_large_bins.yaml -s --antimatter')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_upd_large_bins.yaml -dbshape -s -syst --matter')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_upd_large_bins.yaml -dbshape -s -syst --antimatter')
# os.system('python3 compute_lifetime.py ../Config/2body_analysis_upd_large_bins.yaml -syst --matter')
# os.system('python3 compute_lifetime.py ../Config/2body_analysis_upd_large_bins.yaml -syst --antimatter')



# RUN ANALYSIS WITH SPLITTED MAGNETIC FIELD
# os.system("python3 run_analysis.py ../Config/2body_analysis_B.yaml -t -a -s --matter")
# os.system('python3 compute_lifetime.py ../Config/2body_analysis_B.yaml -dbshape -s -syst --antimatter')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_B.yaml -s -dbshape --matter --syst')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_B.yaml -s -dbshape --antimatter --syst')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_B.yaml -s --matter --syst')
# os.system('python3 compute_blambda.py ../Config/2body_analysis_B.yaml -s --antimatter --syst')