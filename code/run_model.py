import atomica as at
import pandas as pd
from collections import defaultdict
import os
import numpy as np
import constants
import program_instructions

# np.seterr(all='raise')
at.logger.setLevel('DEBUG')
regions = ['SSA', 'SA', 'Other']

coverage_scenarios_max = ['60', '90']

# Folder containing the framework/, databooks/, calibrations/, progset/ and results/ subfolders.
# Update this path to point at your own copy of the analyses folder.
analyses_folder = r'C:\Users\tom.walsh\Projects\A-coupled-model-analyses'
framework_path = os.path.join(analyses_folder, 'framework', 'MNCHdt - stitched framework_20230904.xlsx')

def get_epi_df(scenarios):
    for res in scenarios:
      if "ipmvt: analysis measure" not in res.framework.pars:
        continue
      quantities = res.framework.pars.index[~res.framework.pars["ipmvt: analysis measure"].isna()].tolist()  # Remove parameters where the analysis measure is None or ''
      for quantity in quantities:
        measure = res.framework.pars.at[quantity, "ipmvt: analysis measure"]
        cause = res.framework.pars.at[quantity, "ipmvt: analysis cause"]

        # IPM measures are all aggregated into single year quantities.

        # If the same measure is defined for multiple variables, it will be summed
        for var in res.get_variable(quantity):
          # TODO - This step is slow, would be good to vectorize it somehow e.g. by adding years as columns and unstacking afterwards

          for t, v in zip(np.floor(var.t), var.vals):
            key = (res.name, var.pop.name, int(t), cause, measure)
            if pd.isna(var.timescale):
              records[key] += v
            else:
              records[key] += v / var.timescale * res.dt

          if var.units != at.FrameworkSettings.QUANTITY_TYPE_NUMBER:
            # If quantity is not in number units, average it (within the year)
            for t, count in zip(*np.unique(np.floor(var.t), return_counts=True)):
              key = (res.name, var.pop.name, int(t), cause, measure)
              records[key] /= count

    df = pd.Series(data=records).to_frame(name="Value")
    df.index = df.index.set_names(["Scenario", "Population Name", "Year", "Cause", "Measure"])
    return df

dt = 0.25

# The program coverage TimeSeries and the per-scenario coverage dicts used to build
# at.ProgramInstructions(coverage=...) now live in program_instructions.py.

for region in regions:
  for coverage_scenario_max in coverage_scenarios_max:
    #Programs with region specific values
    coverage_scenarios = pd.read_excel(os.path.join(analyses_folder, 'progset', f'Coverage scenarios_{coverage_scenario_max}%.xlsx'), region)

    coverages = program_instructions.build_coverage_dicts(coverage_scenarios)

    ### Implemented scenarios ###

    #Set paths
    results_dir = os.path.join(analyses_folder, 'results')
    databook_path = os.path.join(analyses_folder, 'databooks', region + '.xlsx')
    calibration_path = os.path.join(analyses_folder, 'calibrations', region + '_calibration.xlsx')

    #import Program and calibrate
    F = at.ProjectFramework(framework_path)
    P = at.Project('MNCHdt', framework=F, do_run=False)
    P.load_databook(databook_path = databook_path, make_default_parset=True, do_run=False)
    P.parsets[0].load_calibration(calibration_path)
    pset = P.load_progbook(os.path.join(analyses_folder, 'progset', 'progset.xlsx'))
    #P = get_latest_project('Optima MNCH', inclusions='MNCHdt', country='Demo')

    # P = fill_assumptions_as_data(P, years=[2020])
    # P.parsets['calibrated'] = P.settings.run_calibration(P)

    #Reconcile program set
    P.settings.update_time_vector(dt=dt)
    P.settings.sim_start = 2015
    P.settings.sim_end = 2050
    # reconciled_pset = run_reconciliation(P, P.parsets[0], P.progsets[0], results_folder=get_desktop_folder(), calib_name = 'Calibrated', recon_name = 'Current spending', reconciliation_plots=False)

    scenarios = ['soc', 'soc_rutf', 'soc_iron_tablets', 'soc_hb_test', 'soc_pe_urine_dipstick', 'soc_ultrasound', 'soc_pe_treatment', 'non_invasive_hb_test',
                 'pe_screening', 'ai_ultrasound', 'soc_scale_up', 'soc_diagnostics_scale_up']

    # scenarios = ['soc']
    scenario_results = []
    for scenario in scenarios:
      #Run scenarios
      print(scenario)
      res = P.run_sim(P.parsets[0], progset=P.progsets[0], progset_instructions=at.ProgramInstructions(start_year=2015, stop_year=2050, coverage = coverages['coverage_' + scenario]), result_name = scenario)
      at.export_results(res, os.path.join(analyses_folder, 'results', scenario + "_programs_" + region + "_" + coverage_scenario_max + ".xlsx"))
      res.export_raw(filename=os.path.join(analyses_folder, 'results', scenario + "_" + region + "_" + coverage_scenario_max + ".xlsx"))
      scenario_results += [res]

    quantities = res.framework.pars.index[~res.framework.pars["ipmvt: analysis measure"].isna()].tolist()
    variables = res.get_variable(quantities[1])
      # d = at.PlotData([res_baseline, res_programs, res_programs_non_diagnosis],['an_severe_anemia_incidence'], ['an_hm'])
      # at.plot_series(d,plot_type='stacked',axis='outputs');
      # plt.title('Anemia treatment');
      # plt.show()




    #Create epi sheet
    records = defaultdict(float)

    epi_df = get_epi_df(scenario_results)
    epi_df = epi_df.query(f"Year >= {constants.analysis_window[0]} & Year <= {constants.analysis_window[1]}")
    epi_df.to_excel(os.path.join(analyses_folder, 'results', "epi_df_" + region + "_" + coverage_scenario_max + ".xlsx"), merge_cells=False)


      # #Create data frame with only outcome parameters
      # df = pd.read_excel('C:\\Users\\Tom.walsh\\Documents\\GitHub\\cat-mnch\\paper\\results\\baseline_results_raw.xlsx')
      # cols = ["Unnamed: 2"]
      #
      # filtered_rows = False
      # for quantity in quantities:
      #     # update filtered_rows with each filter condition
      #     current_filter = (df[cols] == quantity)
      #     filtered_rows |= current_filter
      #
      # filtered_rows = filtered_rows.iloc[:,0]
      # df[filtered_rows]


### Extra scenarios ###
# These exploratory 'scale up' ProgramInstructions objects are not used by the scenario loop above.
# They now live in program_instructions.build_extra_scenario_instructions(coverage_scenarios) -
# call that (with the last-loaded coverage_scenarios DataFrame) if/when they're needed again.



