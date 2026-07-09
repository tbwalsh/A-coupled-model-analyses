import atomica as at
import pandas as pd
from collections import defaultdict
import os
import numpy as np
import constants
from coupled_model import fill_assumptions_as_data


# np.seterr(all='raise')
at.logger.setLevel('DEBUG')
regions = ['SSA']#, 'SSA', 'SA']
coverage_scenarios_max = ['60']

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


####Novel products without displacement####
#Other
ACS_coverage = at.TimeSeries([2022.75, 2023], (np.array([0, 0.9])/dt).tolist())
Az_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
binfantis_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
lung_surfactant_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
neuroprotection_maternal_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
neuroprotection_neonates_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
FCM_pp_iron_moderate_anemia_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
FCM_pp_severe_anemia_coverage = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
c_section_scale_up = at.TimeSeries([2020.75, 2021, 2040], (np.array([0, 0.45, 0.7])/dt).tolist())
hemorrhage_treatment_non_rf_scale_up = at.TimeSeries([2020.75, 2021], (np.array([0, 1])/dt).tolist())
hemorrhage_treatment_rf_scale_up = at.TimeSeries([2020.75, 2021], (np.array([0, 1])/dt).tolist())
ptb_treatment_non_rf_scale_up = at.TimeSeries([2020.75, 2021], (np.array([0, 1])/dt).tolist())
ptb_treatment_rf_scale_up = at.TimeSeries([2020.75, 2021], (np.array([0, 1])/dt).tolist())
referral_facilities_scale_up = at.TimeSeries([2020.75, 2021], (np.array([0, 1])/dt).tolist())
second_trimester_pe_treatment = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())
third_trimester_non_severe_pe_treatment = at.TimeSeries([2020.75, 2021], (np.array([0, 0.9])/dt).tolist())

#Other
mms_scale_up = at.TimeSeries([2015.75, 2016], (np.array([0, 0.9])/dt).tolist())
mms_plus_scale_up = at.TimeSeries([2026, 2027, 2028, 2029, 2030], (np.array([0, 0.1, 0.15, 0.5, 0.8])/dt).tolist())

soc_rutf_scale_up_displacement = at.TimeSeries([2020.75, 2021, 2024, 2025, 2026, 2027], (np.array([0, 0.9, 0.9, 0.6, 0.2, 0])/dt).tolist())
md_rutf_scale_up = at.TimeSeries([2024, 2025, 2026, 2027], (np.array([0, 0.3, 0.7, 0.95])/dt).tolist())

soc_first_trimester_iron_tablets_scale_up_displacement = at.TimeSeries([2020.75, 2021, 2023, 2025], (np.array([0, 0.2, 0.2, 0])/dt).tolist())
soc_second_trimester_iron_tablets_scale_up_displacement = at.TimeSeries([2020.75, 2021, 2023, 2025], (np.array([0, 0.2, 0.2, 0])/dt).tolist())
FCMpregnancy_scale_up = at.TimeSeries([2023, 2025], (np.array([0, 0.9])/dt).tolist())

soc_third_trimester_severe_preclampsia_treatment_scale_up_displacement = at.TimeSeries([2020.75, 2021, 2029, 2032], (np.array([0, 0.5, 0.5, 0])/dt).tolist())
third_trimester_severe_preclampsia_treatment_scale_up = at.TimeSeries([2029, 2032], (np.array([0, 0.95])/dt).tolist())

timeseries = np.arange(2015, 2051, 1).tolist()

for region in regions:
  for coverage_scenario_max in coverage_scenarios_max:
    coverage_scenarios = pd.read_excel(os.path.join(analyses_folder, 'progset', f'Coverage scenarios_{coverage_scenario_max}%.xlsx'), region)

    # SOC - used in SOC scenario
    soc_first_trimester_hb_test_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_rutf_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_rutf_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_iron_tablets_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_iron_tablets_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_iron_tablets_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_iron_tablets_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_third_trimester_severe_preclampsia_treatment_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_severe_preclampsia_treatment_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_ultrasound_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_coverage', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_pe_diagnostic_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_pe_diagnostic_coverage', 2015:2050].div(dt).values.tolist()[0])  # Note diagnostic product but SOC
    soc_third_trimester_pe_diagnostic_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_pe_diagnostic_coverage', 2015:2050].div(dt).values.tolist()[0])

    coverage_soc = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                                                                                           'b.infantisp_routine_nsep_hm': 0, 'lung_surfactant_routine_rds_hm': 0,
                                                                                           'fcm_pp_postpartum_moderate_an_hm': 0,
                                                                                           'fcm_pp_postpartum_severe_an_hm': 0, 'azintrapartum_routine_msep_hm': 0,
                                                                                           'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                                                                           'mms_routine_mnut_hm': mms_scale_up, 'md_rutf_community_19-59_wasting_hm': 0,
                                                                                           'md_rutf_community_6-18_wasting_hm': 0, 'md_rutf_hospitalized_19-59_wasting_hm': 0,
                                                                                           'md_rutf_hospitalized_6-18_wasting_hm': 0, 'azpregnancy_routine_msep_hm': 0,
                                                                                           'proxy_preterm_treatment___non_referral_routine_pt_hm': 0,
                                                                                           'mms__routine_mnut_hm': 0,
                                                                                           'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                                                                           'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                                                                           'preeclampsia_novel_therapeutic_drug_second_trimester_preecl_hm': 0,
                                                                                           'preeclampsia_novel_therapeutic_drug_third_trimester_non_severe_preecl_hm': 0,
                                                                                           'preeclampsia_novel_therapeutic_drug_third_trimester_severe_preecl_hm': 0,
                                                                                           'b.infantistx_routine_wasting_hm': 0, 'acs_routine_rds_hm':0,
                                                                                           'proxy_hemorrhage_treatment___referral_routine_hm_hm': 0,
                                                                                           'ifa_routine_mnut_hm': 0,
                                                                                           'neuroprotection_candidates_routine_maternal_ol_hm': 0,
                                                                                           'neuroprotection_candidates_routine_neonates_asx_hm': 0,
                                                                                           'proxy_hemorrhage_treatment___non_referral_routine_hm_hm': 0,
                                                                                           'fcmpregnancy_moderate_second_trimester_an_hm': 0,
                                                                                           'fcmpregnancy_severe_second_trimester_an_hm': 0,
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': 0,
                                                                                           'proxy_preterm_treatment___referral_routine_pt_hm': 0,
                                                                                           'proxy_referral_facility_referral_facility_anc_hm': 0,
                                                                                           'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                                                                                           'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': 0
                                                                                           }


    coverages = {'coverage_soc':coverage_soc}

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

    P = fill_assumptions_as_data(P, years=[2020])
    # P.parsets['calibrated'] = P.settings.run_calibration(P)

    #Reconcile program set
    P.settings.update_time_vector(dt=dt)
    P.settings.sim_start = 2015
    P.settings.sim_end = 2050
    #reconciled_pset = run_reconciliation(P, P.parsets[0], P.progsets[0], results_folder=get_desktop_folder(), calib_name = 'Calibrated', recon_name = 'Current spending', reconciliation_plots=False)

    scenarios = ['soc']

    # scenarios = ['soc']
    scenario_results = []
    for scenario in scenarios:
      #Run scenarios
      print(scenario)
      res = P.run_sim(P.parsets[0], progset=P.progsets[0], progset_instructions=at.ProgramInstructions(start_year=2015, stop_year=2050, coverage = coverages['coverage_' + scenario]), result_name = scenario)
      res.export_raw(filename=os.path.join(analyses_folder, 'results', scenario + "_" + region + ".xlsx"))
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

if region == 'SA':
    region = 'South Asia'
    acronym = 'SA'
elif region == 'SSA':
    region = 'Sub-Saharan Africa'
    acronym = 'SSA'
else:
    region = 'Other LMICs'
    acronym = 'Other'
df = pd.read_excel(os.path.join(analyses_folder, 'calibrations', 'Calibration_vals_2023.xlsx'), index_col = 'code_name', sheet_name = acronym)
model_data = pd.read_excel(os.path.join(analyses_folder, 'results', f'soc_{acronym}.xlsx'), index_col = 2)

an_bleeding_related_anemia = df._get_value('an_bleeding_related_anemia','Middle')
pt_prop_sga_pt = df._get_value('pt_incident_cases_sga_pt','Middle')
pt_prop_sga_term = df._get_value('pt_incident_cases_sga_term','Middle')
pt_incident_cases_pt=df._get_value('pt_incident_cases_pt','Middle')
pt_mortalities_pt=df[(df['Region']==region)]._get_value('pt_mortalities_pt','Middle')
pt_incident_cases_sga_pt=df[(df['Region']==region)]._get_value('pt_incident_cases_sga_pt','Middle')
pt_incident_cases_sga_term=df[(df['Region']==region)]._get_value('pt_incident_cases_sga_term','Middle')
hm_incident_cases=df[(df['Region']==region)]._get_value('hm_incident_cases','Middle')
hm_mortalities=df[(df['Region']==region)]._get_value('hm_mortalities','Middle')
ol_incident_cases=df[(df['Region']==region)]._get_value('ol_incident_cases','Middle')
ol_mortalities=df[(df['Region']==region)]._get_value('ol_mortalities','Middle')
msep_incidence=df[(df['Region']==region)]._get_value('msep_incidence','Middle')
msep_mortalities=df[(df['Region']==region)]._get_value('msep_mortalities','Middle')
asx_num_live_births=df[(df['Region']==region)]._get_value('asx_num_live_births','Middle')
asx_stillbirth_incidence=df[(df['Region']==region)]._get_value('asx_stillbirth_incidence','Middle')
rds_mortalities=df[(df['Region']==region)]._get_value('rds_mortalities','Middle')
asx_birth_Asphyxia_incidence=df[(df['Region']==region)]._get_value('asx_birth_Asphyxia_incidence','Middle')
asx_encephalopathy_mortality=df[(df['Region']==region)]._get_value('asx_encephalopathy_mortality','Middle')
nsep_incident_cases=df[(df['Region']==region)]._get_value('nsep_incident_cases','Middle')
nsep_mortalities=df[(df['Region']==region)]._get_value('nsep_mortalities','Middle')
oth_inf_incident_cases = df[(df['Region']==region)]._get_value('oth_inf_incident_cases','Middle')
oth_inf_mortalities = df[(df['Region']==region)]._get_value('oth_inf_mortalities','Middle')
preecl_incident_cases = df[(df['Region']==region)]._get_value('preecl_incident_cases','Middle')
preecl_mortalities = df[(df['Region']==region)]._get_value('preecl_mortalities','Middle')
wasting_dalys = df[(df['Region']==region)]._get_value('wasting_dalys','Middle')
stunting_dalys = df[(df['Region']==region)]._get_value('stunting_dalys','Middle')
wasting_prop_1 = df[(df['Region']==region)]._get_value('wasting_prop_wasting_occurs_1','Middle')


mult_factor = []
mult_factor += [an_bleeding_related_anemia/model_data.loc['an_prop_healthy_bleeding_related_anemia_pp', 2019]]
mult_factor += [pt_incident_cases_sga_pt/model_data.loc['pt_prop_sga_preterm', 2019]]
mult_factor += [pt_incident_cases_sga_term/model_data.loc['pt_prop_sga_term', 2019]]
mult_factor += [pt_incident_cases_pt/(model_data.loc['pt_incident_cases_pt', 2019]+model_data.loc['pt_incident_cases_pt_sga', 2019])]
mult_factor += [pt_mortalities_pt/model_data.loc['pt_mortalities_pt', 2019]]
mult_factor += [hm_incident_cases/model_data.loc['hm_incident_cases', 2019]]
mult_factor += [hm_mortalities/model_data.loc['hm_mortalities', 2019]]
mult_factor += [ol_incident_cases/model_data.loc['ol_incident_cases', 2019]]
mult_factor += [ol_mortalities/model_data.loc['ol_mortalities', 2019]]
mult_factor += [msep_incidence/model_data.loc['msep_incidence', 2019]]
mult_factor += [msep_mortalities/model_data.loc['msep_mortalities', 2019]]
mult_factor += [preecl_incident_cases/model_data.loc['preecl_incident_cases', 2019]]
mult_factor += [preecl_mortalities/model_data.loc['preecl_mortalities', 2019]]
mult_factor += [asx_num_live_births/model_data.loc['asx_num_live_births', 2019]]
mult_factor += [asx_stillbirth_incidence/model_data.loc['asx_stillbirth_incidence', 2019]]
mult_factor += [rds_mortalities/model_data.loc['rds_mortalities', 2019]]
mult_factor += [asx_birth_Asphyxia_incidence/model_data.loc['asx_birth_Asphyxia_incidence', 2019]]
mult_factor += [asx_encephalopathy_mortality/model_data.loc['asx_encephalopathy_mortality', 2019]]
mult_factor += [nsep_incident_cases/model_data.loc['nsep_incident_cases', 2019]]
mult_factor += [nsep_mortalities/model_data.loc['nsep_mortalities', 2019]]
#mult_factor += [oth_inf_mortalities/name.model.get_pop('oth_inf_hm').get_par('oth_inf_mortalities').vals[20]]
mult_factor += [wasting_prop_1/model_data.loc['wasting_prop_wasting_occurs_1', 2019]] #calibration for first timestep
mult_factor += [wasting_dalys/model_data.loc['wasting_dalys', 2019]]
mult_factor += [stunting_dalys/model_data.loc['stunting_dalys', 2019]]
mult_factor_dict = {
                    'pt_prop_sga_preterm_no_risk': mult_factor[1],
                    'pt_prop_sga_term_no_risk': mult_factor[2],
                    'pt_incident_cases_pt': mult_factor[3],
                    'hm_incident_cases': mult_factor[5],
                    'ol_incident_cases': mult_factor[7],
                    'msep_incidence': mult_factor[9],
                    'preecl_incident_cases': mult_factor[11],
                    'pt_mortalities_pt': mult_factor [4],
                    'hm_mortalities': mult_factor [6],
                    'ol_mortalities': mult_factor[8],
                    'msep_mortalities': mult_factor [10],
                    'preecl_mortalities': mult_factor [12],
                    'asx_num_live_births': mult_factor [13],
                    'an_prop_healthy_bleeding_related_anemia_pp': mult_factor[0],
                    'asx_stillbirth_incidence': mult_factor [14],
                    'rds_mortalities': mult_factor [15],
                    'asx_birth_Asphyxia_incidence': mult_factor [16],
                    'nsep_incident_cases': mult_factor[18],
                    'wasting_incidence_1': mult_factor[20],
                    'stunting_dalys': mult_factor[22],
                    'asx_encephalopathy_mortality': mult_factor [17],
                    'nsep_mortalities': mult_factor [19],
                    'wasting_dalys': mult_factor [21]}

print(mult_factor_dict)