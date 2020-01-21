# a test of using matlab functions from sysplotter in Python
import matlab.engine
eng = matlab.engine.start_matlab()

home = "/home/qkonyn/Programming_Projects/GeometricSystemPlotter/"
eng.addpath(home + "UserFiles/GenericUser/Systems/")
eng.addpath(home + "ProgramFiles/Physics/LowReynoldsRFT/Discrete_links/")

s = eng.sysf_two_link_lowRe()

shape = 2
backwards = [False, False]
eng.LowRE_connection_discrete(s["geometry"], s["physics"], shape, backwards)


eng.quit()
