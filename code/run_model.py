import atomica as at
from at_tools import get_latest_project, verify_model, get_desktop_folder, get_apps_folder, get_paths, stitch_frameworks, fill_assumptions_as_data, run_reconciliation
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import gcat as gc
from os import sep
import numpy as np
import sciris as sc
import mnch

# np.seterr(all='raise')
at.logger.setLevel('DEBUG')
regions = ['Afghanistan']

coverage_scenarios_max = ['60', '90']


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


#Programs without region specific values
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
    #Programs with region specific values
    coverage_scenarios = pd.read_excel(f"C:\\Users\\Tom.walsh\\OneDrive - Burnet Institute\\Desktop\\Tom Walsh PHD work\\PhD work\\Papers\\Paper 3 PPH\\Coverage_scenarios.xlsx", region)

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

    # SOC scale up - used in SOC scale up scenario
    soc_first_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_rutf_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_rutf_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_iron_tablets_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_iron_tablets_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_iron_tablets_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_iron_tablets_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_third_trimester_severe_preclampsia_treatment_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_severe_preclampsia_treatment_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_ultrasound_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_pe_diagnostic_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_pe_diagnostic_scale_up', 2015:2050].div(dt).values.tolist()[0])  # Note diagnostic product but SOC
    soc_third_trimester_pe_diagnostic_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_pe_diagnostic_scale_up', 2015:2050].div(dt).values.tolist()[0])  # Note diagnostic product but SOC

    #### Diagnostics ####
    ####Novel products without displacement####
    preeclampsia_screening_prevention_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'preeclampsia_screening_prevention_coverage', 2015:2050].div(dt).values.tolist()[0])

    # Displacement
    soc_ultrasound_coverage_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_coverage_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_ultrasound_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    ai_ultrasound_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'ai_ultrasound_scale_up', 2015:2050].div(dt).values.tolist()[0])  # Note diagnostic product (90% in actual scenario)

    # redo these
    soc_first_trimester_hb_test_coverage_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_coverage_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_hb_test_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    novel_first_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'novel_first_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])  # Note diagnostic product
    novel_second_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'novel_second_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0]) # Note diagnostic product

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

    coverage_soc_rutf = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                                                                                           'b.infantisp_routine_nsep_hm': 0, 'lung_surfactant_routine_rds_hm': 0,
                                                                                           'fcm_pp_postpartum_moderate_an_hm': 0, 'fcm_pp_postpartum_severe_an_hm': 0, 'azintrapartum_routine_msep_hm': 0, 'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                                                                           'mms_routine_mnut_hm': mms_scale_up, 'md_rutf_community_19-59_wasting_hm': 0, 'md_rutf_community_6-18_wasting_hm': 0, 'md_rutf_hospitalized_19-59_wasting_hm': 0,
                                                                                           'md_rutf_hospitalized_6-18_wasting_hm': 0, 'azpregnancy_routine_msep_hm': 0, 'proxy_preterm_treatment___non_referral_routine_pt_hm': 0,
                                                                                           'mms__routine_mnut_hm': 0, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                                                                           'preeclampsia_novel_therapeutic_drug_second_trimester_preecl_hm': 0, 'preeclampsia_novel_therapeutic_drug_third_trimester_non_severe_preecl_hm': 0,
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

    coverage_soc_iron_tablets={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                                                                                           'b.infantisp_routine_nsep_hm': 0, 'lung_surfactant_routine_rds_hm': 0,
                                                                                           'fcm_pp_postpartum_moderate_an_hm': 0, 'fcm_pp_postpartum_severe_an_hm': 0, 'azintrapartum_routine_msep_hm': 0, 'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                                                                           'mms_routine_mnut_hm': mms_scale_up, 'md_rutf_community_19-59_wasting_hm': 0, 'md_rutf_community_6-18_wasting_hm': 0, 'md_rutf_hospitalized_19-59_wasting_hm': 0,
                                                                                           'md_rutf_hospitalized_6-18_wasting_hm': 0, 'azpregnancy_routine_msep_hm': 0, 'proxy_preterm_treatment___non_referral_routine_pt_hm': 0,
                                                                                           'mms__routine_mnut_hm': 0, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                                                                           'preeclampsia_novel_therapeutic_drug_second_trimester_preecl_hm': 0, 'preeclampsia_novel_therapeutic_drug_third_trimester_non_severe_preecl_hm': 0,
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

    coverage_soc_hb_test = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
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

    coverage_soc_pe_urine_dipstick={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
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

    coverage_soc_ultrasound ={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
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

    coverage_soc_pe_treatment = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
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

    coverage_non_invasive_hb_test={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up_displacement, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up_displacement,
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
                                                                                           'noninvasive_hb_diagnostic_first_trimester_an_hm': novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': 0
                                                                                           }

    coverage_pe_screening = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up_displacement, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
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
                                                                                           'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': third_trimester_severe_preclampsia_treatment_scale_up,
                                                                                           'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': 0
                                                                                           }

    coverage_ai_ultrasound={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up_displacement, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
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
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                                                                                           }

    coverage_soc_scale_up={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
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

    coverage_soc_diagnostics_scale_up={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                                                                                           'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                                                                                           'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                                                                           'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up_displacement, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up_displacement, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up_displacement,
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
                                                                                           'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':preeclampsia_screening_prevention_coverage,
                                                                                           'noninvasive_hb_diagnostic_first_trimester_an_hm': novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                                                                           'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                                                                                           }

    coverages = {'coverage_soc':coverage_soc, 'coverage_soc_rutf':coverage_soc_rutf, 'coverage_soc_iron_tablets':coverage_soc_iron_tablets,
                 'coverage_soc_hb_test':coverage_soc_hb_test, 'coverage_soc_pe_urine_dipstick':coverage_soc_pe_urine_dipstick, 'coverage_soc_ultrasound':coverage_soc_ultrasound, 'coverage_soc_pe_treatment':coverage_soc_pe_treatment,
                 'coverage_non_invasive_hb_test':coverage_non_invasive_hb_test, 'coverage_pe_screening':coverage_pe_screening, 'coverage_ai_ultrasound':coverage_ai_ultrasound,
                 'coverage_soc_scale_up':coverage_soc_scale_up, 'coverage_soc_diagnostics_scale_up':coverage_soc_diagnostics_scale_up}

    ### Implemented scenarios ###

    #Set paths
    results_dir = 'C:\\Users\\Tom.walsh\\OneDrive - Burnet Institute\\Desktop\\Tom Walsh PHD work\\PhD work\\Papers\\Paper 3 PPH\\Atomica books\\'
    databook_path = 'C:\\Users\\Tom.walsh\\OneDrive - Burnet Institute\\Desktop\\Tom Walsh PHD work\\PhD work\\Papers\\Paper 3 PPH\\Atomica books\\Country databooks\\' + region + '.xlsx'
    calibration_path = 'C:\\Users\\Tom.walsh\\OneDrive - Burnet Institute\\Desktop\\Tom Walsh PHD work\\PhD work\\Papers\\Paper 3 PPH\\Atomica books\\' + region + '_calibration.xlsx'

    #import Program and calibrate
    F = at.ProjectFramework(f'C:\\Users\\Tom.walsh\\OneDrive - Burnet Institute\\Desktop\\Tom Walsh PHD work\\PhD work\\Papers\\Paper 3 PPH\\Atomica books\\MNCHdt - stitched framework_20250819.xlsx')
    P = at.Project('MNCHdt', framework=F, do_run=False)
    P.load_databook(databook_path = databook_path, make_default_parset=True, do_run=False)
    P.parsets[0].load_calibration(calibration_path)
    pset = P.load_progbook('C:\\Users\\Tom.walsh\\Documents\\GitHub\\cat-mnch\\paper\\progset\\progset.xlsx')
    #P = get_latest_project('Optima MNCH', inclusions='MNCHdt', country='Demo')

    P = fill_assumptions_as_data(P, years=[2020])
    # P.parsets['calibrated'] = P.settings.run_calibration(P)

    #Reconcile program set
    P.settings.update_time_vector(dt=dt)
    P.settings.sim_start = 2015
    P.settings.sim_end = 2050
    reconciled_pset = run_reconciliation(P, P.parsets[0], P.progsets[0], results_folder=get_desktop_folder(), calib_name = 'Calibrated', recon_name = 'Current spending', reconciliation_plots=False)

    scenarios = ['soc', 'soc_rutf', 'soc_iron_tablets', 'soc_hb_test', 'soc_pe_urine_dipstick', 'soc_ultrasound', 'soc_pe_treatment', 'non_invasive_hb_test',
                 'pe_screening', 'ai_ultrasound', 'soc_scale_up', 'soc_diagnostics_scale_up']

    # scenarios = ['soc']
    scenario_results = []
    for scenario in scenarios:
      #Run scenarios
      print(scenario)
      res = P.run_sim(P.parsets[0], progset=P.progsets[0], progset_instructions=at.ProgramInstructions(start_year=2015, stop_year=2050, coverage = coverages['coverage_' + scenario]), result_name = scenario)
      at.export_results(res, "C:\\Users\\Tom.walsh\\Documents\\GitHub\\cat-mnch\\paper\\results\\" + scenario + "_programs_" + region + "_" + coverage_scenario_max + ".xlsx")
      res.export_raw(filename="C:\\Users\\Tom.walsh\\Documents\\GitHub\\cat-mnch\\paper\\results\\" + scenario + "_" + region + "_" + coverage_scenario_max + ".xlsx")
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
    epi_df = epi_df.query(f"Year >= {mnch.analysis_window[0]} & Year <= {mnch.analysis_window[1]}")
    epi_df.to_excel("C:\\Users\\Tom.walsh\\Documents\\GitHub\\cat-mnch\\paper\\results\\epi_df_" + region + "_" + coverage_scenario_max + ".xlsx", merge_cells=False)


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
instructions_scaleup = at.ProgramInstructions(start_year=2020, stop_year = 2040, coverage={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up_displacement,
                                                                                       'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up_displacement,
                                                                                       'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up_displacement, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                                                                       'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up_displacement, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up_displacement, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up_displacement, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up_displacement,
                                                                                           'b.infantisp_routine_nsep_hm':binfantis_coverage, 'lung_surfactant_routine_rds_hm':lung_surfactant_coverage, 'fcm_pp_postpartum_moderate_an_hm':FCM_pp_iron_moderate_anemia_coverage, 'fcm_pp_postpartum_severe_an_hm':FCM_pp_severe_anemia_coverage, 'azintrapartum_routine_msep_hm':Az_coverage, 'c_section___obstructed_labor_routine_ol_hm':c_section_scale_up,
                                                                                           'mms_routine_mnut_hm':mms_scale_up, 'md_rutf_community_19-59_wasting_hm':md_rutf_scale_up, 'md_rutf_community_6-18_wasting_hm':md_rutf_scale_up, 'md_rutf_hospitalized_19-59_wasting_hm':md_rutf_scale_up, 'md_rutf_hospitalized_6-18_wasting_hm':md_rutf_scale_up, 'azpregnancy_routine_msep_hm':Az_coverage,
                                                                                           'proxy_preterm_treatment___non_referral_routine_pt_hm':ptb_treatment_non_rf_scale_up, 'mms__routine_mnut_hm':mms_plus_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm':novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm':novel_second_trimester_hb_test_scale_up, 'preeclampsia_novel_therapeutic_drug_second_trimester_preecl_hm':second_trimester_pe_treatment,
                                                                                           'preeclampsia_novel_therapeutic_drug_third_trimester_non_severe_preecl_hm':third_trimester_non_severe_pe_treatment, 'preeclampsia_novel_therapeutic_drug_third_trimester_severe_preecl_hm':third_trimester_severe_preclampsia_treatment_scale_up, 'b.infantistx_routine_wasting_hm':binfantis_coverage, 'acs_routine_rds_hm':ACS_coverage, 'proxy_hemorrhage_treatment___referral_routine_hm_hm':hemorrhage_treatment_rf_scale_up,
                                                                                           'ifa_routine_mnut_hm':mms_plus_scale_up, 'neuroprotection_candidates_routine_maternal_ol_hm':neuroprotection_neonates_coverage, 'neuroprotection_candidates_routine_neonates_asx_hm':neuroprotection_maternal_coverage, 'proxy_hemorrhage_treatment___non_referral_routine_hm_hm':hemorrhage_treatment_non_rf_scale_up, 'fcmpregnancy_moderate_second_trimester_an_hm':FCMpregnancy_scale_up,
                                                                                           'fcmpregnancy_severe_second_trimester_an_hm':FCMpregnancy_scale_up, 'ai_enabled_ultrasound_routine_anc_hm':ai_ultrasound_scale_up, 'proxy_preterm_treatment___referral_routine_pt_hm':ptb_treatment_rf_scale_up, 'proxy_referral_facility_referral_facility_anc_hm':referral_facilities_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm':novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                                                                           'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':preeclampsia_screening_prevention_coverage})

instructions_scaleup_non_diagnosis = at.ProgramInstructions(start_year=2020, stop_year = 2040, coverage={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up_displacement, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up_displacement,
                                                                                       'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up_displacement,
                                                                                       'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up_displacement, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                                                                       'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up_displacement, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
                                                                                           'b.infantisp_routine_nsep_hm':binfantis_coverage, 'lung_surfactant_routine_rds_hm':lung_surfactant_coverage, 'fcm_pp_postpartum_moderate_an_hm':FCM_pp_iron_moderate_anemia_coverage, 'fcm_pp_postpartum_severe_an_hm':FCM_pp_severe_anemia_coverage, 'azintrapartum_routine_msep_hm':Az_coverage, 'c_section___obstructed_labor_routine_ol_hm':c_section_scale_up,
                                                                                           'mms_routine_mnut_hm':mms_scale_up, 'md_rutf_community_19-59_wasting_hm':md_rutf_scale_up, 'md_rutf_community_6-18_wasting_hm':md_rutf_scale_up, 'md_rutf_hospitalized_19-59_wasting_hm':md_rutf_scale_up, 'md_rutf_hospitalized_6-18_wasting_hm':md_rutf_scale_up, 'azpregnancy_routine_msep_hm':Az_coverage,
                                                                                           'proxy_preterm_treatment___non_referral_routine_pt_hm':ptb_treatment_non_rf_scale_up, 'mms__routine_mnut_hm':mms_plus_scale_up, 'preeclampsia_novel_therapeutic_drug_second_trimester_preecl_hm':second_trimester_pe_treatment,
                                                                                           'preeclampsia_novel_therapeutic_drug_third_trimester_non_severe_preecl_hm':third_trimester_non_severe_pe_treatment, 'preeclampsia_novel_therapeutic_drug_third_trimester_severe_preecl_hm':third_trimester_severe_preclampsia_treatment_scale_up, 'b.infantistx_routine_wasting_hm':binfantis_coverage, 'acs_routine_rds_hm':ACS_coverage, 'proxy_hemorrhage_treatment___referral_routine_hm_hm':hemorrhage_treatment_rf_scale_up,
                                                                                           'ifa_routine_mnut_hm':mms_plus_scale_up, 'neuroprotection_candidates_routine_maternal_ol_hm':neuroprotection_neonates_coverage, 'neuroprotection_candidates_routine_neonates_asx_hm':neuroprotection_maternal_coverage, 'proxy_hemorrhage_treatment___non_referral_routine_hm_hm':hemorrhage_treatment_non_rf_scale_up, 'fcmpregnancy_moderate_second_trimester_an_hm':FCMpregnancy_scale_up,
                                                                                           'fcmpregnancy_severe_second_trimester_an_hm':FCMpregnancy_scale_up, 'proxy_preterm_treatment___referral_routine_pt_hm':ptb_treatment_rf_scale_up, 'proxy_referral_facility_referral_facility_anc_hm':referral_facilities_scale_up,
                                                                                                         #Set diagnostics to 0
                                                                                          'noninvasive_hb_diagnostic_first_trimester_an_hm':novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up, 'ai_enabled_ultrasound_routine_anc_hm':0, 'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':0})


instructions_scaleup_diagnosis_only = at.ProgramInstructions(start_year=2020, stop_year = 2040, coverage={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                                                                                       'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                                                                                       'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                                                                       'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
                                                                                       'b.infantisp_routine_nsep_hm': 0, 'lung_surfactant_routine_rds_hm': 0,
                                                                                       'fcm_pp_postpartum_moderate_an_hm': 0,
                                                                                       'fcm_pp_postpartum_severe_an_hm': 0, 'azintrapartum_routine_msep_hm': 0,
                                                                                       'c_section___obstructed_labor_routine_ol_hm': 0,
                                                                                       'mms_routine_mnut_hm': 0, 'md_rutf_community_19-59_wasting_hm': 0,
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
                                                                                       'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':preeclampsia_screening_prevention_coverage,
                                                                                       'noninvasive_hb_diagnostic_first_trimester_an_hm':novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                                                                       'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                                                                                       })



