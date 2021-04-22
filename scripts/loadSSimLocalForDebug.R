e$PackageDirectory = "C:/Users/HughesJo/Documents/SyncroSim/Packages/ROFSim"
e$ProjectId = 1
e$ScenarioId = 30
myLib <- ssimLibrary("C:/Users/HughesJo/Documents/InitialWork/OntarioFarNorth/RoFModel/UI/ROFSim")
mySce <- scenario(myLib,e$ScenarioId)

GLOBAL_Library = myLib#ssimLibrary(session = GLOBAL_Session)
GLOBAL_Project = project(GLOBAL_Library, project = as.integer(e$ProjectId))
GLOBAL_Scenario = scenario(GLOBAL_Library, scenario = as.integer(e$ScenarioId))
GLOBAL_RunControl = GetDataSheetExpectData("ROFSim_RunControl", GLOBAL_Scenario)
GLOBAL_MaxIteration = GetSingleValueExpectData(GLOBAL_RunControl, "MaximumIteration")
GLOBAL_MinIteration = GetSingleValueExpectData(GLOBAL_RunControl, "MinimumIteration")
GLOBAL_MinTimestep = GetSingleValueExpectData(GLOBAL_RunControl, "MinimumTimestep")
GLOBAL_MaxTimestep = GetSingleValueExpectData(GLOBAL_RunControl, "MaximumTimestep")
GLOBAL_TotalIterations = (GLOBAL_MaxIteration - GLOBAL_MinIteration + 1)
GLOBAL_TotalTimesteps = (GLOBAL_MaxTimestep - GLOBAL_MinTimestep + 1)

#use local dev version of caribou metrics for debugging.
detach("package:caribouMetrics", unload=TRUE)
devtools::document();devtools::load_all()
