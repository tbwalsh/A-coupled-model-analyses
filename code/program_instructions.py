import atomica as at
import numpy as np

# This module builds the atomica program coverage dictionaries / ProgramInstructions used by run_model.py.
# It was previously defined inline in run_model.py; moved here so run_model.py only has to call it.

dt = 0.25

timeseries = np.arange(2015, 2051, 1).tolist()

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


def build_coverage_dicts(coverage_scenarios):
    """
    Build the {'coverage_<scenario>': {program: TimeSeries/value}} lookup used to construct
    each scenario's at.ProgramInstructions(coverage=...) for one region/coverage_scenario_max.

    :param coverage_scenarios: DataFrame loaded from 'Coverage scenarios_<60/90>%.xlsx' for one region
    :return: dict of {'coverage_soc': {...}, 'coverage_soc_rutf': {...}, ...}
    """

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
                     'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                     'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                     'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                     'ai_enabled_ultrasound_routine_anc_hm': 0,
                     'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                     'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                     'ai_enabled_ultrasound_routine_anc_hm': 0
                     }

    coverage_soc_rutf = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                          'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                          'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                          'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                          'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                          'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                          'ai_enabled_ultrasound_routine_anc_hm': 0,
                          'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                          'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                          'ai_enabled_ultrasound_routine_anc_hm': 0
                          }

    coverage_soc_iron_tablets={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                               'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                               'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                               'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                               'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                               'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                               'ai_enabled_ultrasound_routine_anc_hm': 0,
                               'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                               'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                               'ai_enabled_ultrasound_routine_anc_hm': 0
                               }

    coverage_soc_hb_test = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                             'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                             'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                             'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
                             'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                             'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                             'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                             'ai_enabled_ultrasound_routine_anc_hm': 0,
                             'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                             'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                             'ai_enabled_ultrasound_routine_anc_hm': 0
                             }

    coverage_soc_pe_urine_dipstick={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                     'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                     'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                     'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                                     'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                     'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                     'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                     'ai_enabled_ultrasound_routine_anc_hm': 0,
                                     'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                                     'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                                     'ai_enabled_ultrasound_routine_anc_hm': 0
                                     }

    coverage_soc_ultrasound ={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                               'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                               'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                               'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                               'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                               'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                               'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                               'ai_enabled_ultrasound_routine_anc_hm': 0,
                               'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                               'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                               'ai_enabled_ultrasound_routine_anc_hm': 0
                               }

    coverage_soc_pe_treatment = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                  'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                  'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                  'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                                  'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                  'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                  'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                  'ai_enabled_ultrasound_routine_anc_hm': 0,
                                  'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                                  'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                                  'ai_enabled_ultrasound_routine_anc_hm': 0
                                  }

    coverage_non_invasive_hb_test={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                                    'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                                    'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                                    'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up_displacement, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up_displacement,
                                    'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                    'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                    'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                    'ai_enabled_ultrasound_routine_anc_hm': 0,
                                    'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                                    'noninvasive_hb_diagnostic_first_trimester_an_hm': novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                    'ai_enabled_ultrasound_routine_anc_hm': 0
                                    }

    coverage_pe_screening = {'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                              'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                              'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                              'soc_ultrasound_routine_anc_hm':soc_ultrasound_coverage, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up_displacement, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                              'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                              'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                              'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                              'ai_enabled_ultrasound_routine_anc_hm': 0,
                              'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': third_trimester_severe_preclampsia_treatment_scale_up,
                              'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                              'ai_enabled_ultrasound_routine_anc_hm': 0
                              }

    coverage_ai_ultrasound={'soc_rutf_community_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_coverage, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_coverage,
                             'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage,
                             'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_coverage, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_coverage, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_coverage, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_coverage,
                             'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up_displacement, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_coverage, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_coverage, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_coverage,
                             'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                             'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                             'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                             'ai_enabled_ultrasound_routine_anc_hm': 0,
                             'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                             'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                             'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                             }

    coverage_soc_scale_up={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                            'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                            'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                            'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up,
                            'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                            'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                            'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                            'ai_enabled_ultrasound_routine_anc_hm': 0,
                            'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm': 0,
                            'noninvasive_hb_diagnostic_first_trimester_an_hm': 0, 'noninvasive_hb_diagnostic_second_trimester_an_hm': 0,
                            'ai_enabled_ultrasound_routine_anc_hm': 0
                            }

    coverage_soc_diagnostics_scale_up={'soc_rutf_community_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_community_6-18_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_19-59_wasting_hm':soc_rutf_scale_up, 'soc_rutf_hospitalized_6-18_wasting_hm':soc_rutf_scale_up,
                                        'soc_iron_tablets_mild_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_mild_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_moderate_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up,
                                        'soc_iron_tablets_severe_first_trimester_an_hm':soc_first_trimester_iron_tablets_scale_up, 'soc_iron_tablets_severe_second_trimester_an_hm':soc_second_trimester_iron_tablets_scale_up, 'soc_preeclampsia_urine_dipstick_test_second_trimester_preecl_hm':soc_second_trimester_pe_diagnostic_scale_up, 'soc_preeclampsia_urine_dipstick_test_third_trimester_preecl_hm':soc_third_trimester_pe_diagnostic_scale_up,
                                        'soc_ultrasound_routine_anc_hm':soc_ultrasound_scale_up_displacement, 'soc_preeclampsia_treatment_third_trimester_severe_preecl_hm':soc_third_trimester_severe_preclampsia_treatment_scale_up, 'soc_hb_test_first_trimester_an_hm':soc_first_trimester_hb_test_scale_up_displacement, 'soc_hb_test_second_trimester_an_hm':soc_second_trimester_hb_test_scale_up_displacement,
                                        'c_section___obstructed_labor_routine_ol_hm': c_section_scale_up,
                                        'mms_routine_mnut_hm': mms_scale_up, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                        'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                        'ai_enabled_ultrasound_routine_anc_hm': 0,
                                        'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':preeclampsia_screening_prevention_coverage,
                                        'noninvasive_hb_diagnostic_first_trimester_an_hm': novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                        'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                                        }

    coverages = {'coverage_soc':coverage_soc, 'coverage_soc_rutf':coverage_soc_rutf, 'coverage_soc_iron_tablets':coverage_soc_iron_tablets,
                 'coverage_soc_hb_test':coverage_soc_hb_test, 'coverage_soc_pe_urine_dipstick':coverage_soc_pe_urine_dipstick, 'coverage_soc_ultrasound':coverage_soc_ultrasound, 'coverage_soc_pe_treatment':coverage_soc_pe_treatment,
                 'coverage_non_invasive_hb_test':coverage_non_invasive_hb_test, 'coverage_pe_screening':coverage_pe_screening, 'coverage_ai_ultrasound':coverage_ai_ultrasound,
                 'coverage_soc_scale_up':coverage_soc_scale_up, 'coverage_soc_diagnostics_scale_up':coverage_soc_diagnostics_scale_up}

    return coverages


def build_extra_scenario_instructions(coverage_scenarios):
    """
    Build the three exploratory 'scale up' ProgramInstructions objects ('Extra scenarios' in the
    original run_model.py). These are not currently called by run_model.py's main scenario loop -
    kept here for reference/future use.

    :param coverage_scenarios: DataFrame loaded from 'Coverage scenarios_<60/90>%.xlsx' for one region
    :return: dict of {'instructions_scaleup': ..., 'instructions_scaleup_non_diagnosis': ..., 'instructions_scaleup_diagnosis_only': ...}
    """

    soc_rutf_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_rutf_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_iron_tablets_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_iron_tablets_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_iron_tablets_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_iron_tablets_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_pe_diagnostic_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_pe_diagnostic_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_third_trimester_pe_diagnostic_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_pe_diagnostic_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_ultrasound_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_ultrasound_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_ultrasound_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    soc_first_trimester_hb_test_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_first_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_second_trimester_hb_test_scale_up_displacement = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_second_trimester_hb_test_scale_up_displacement', 2015:2050].div(dt).values.tolist()[0])
    soc_third_trimester_severe_preclampsia_treatment_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'soc_third_trimester_severe_preclampsia_treatment_scale_up', 2015:2050].div(dt).values.tolist()[0])
    novel_first_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'novel_first_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    novel_second_trimester_hb_test_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'novel_second_trimester_hb_test_scale_up', 2015:2050].div(dt).values.tolist()[0])
    preeclampsia_screening_prevention_coverage = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'preeclampsia_screening_prevention_coverage', 2015:2050].div(dt).values.tolist()[0])
    ai_ultrasound_scale_up = at.TimeSeries(timeseries, coverage_scenarios.loc[coverage_scenarios["Product"] == 'ai_ultrasound_scale_up', 2015:2050].div(dt).values.tolist()[0])

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
                                                                                                       'c_section___obstructed_labor_routine_ol_hm': 0,
                                                                                                       'mms_routine_mnut_hm': 0, 'noninvasive_hb_diagnostic_first_trimester_an_hm': 0,
                                                                                                       'noninvasive_hb_diagnostic_second_trimester_an_hm':0,
                                                                                                       'ai_enabled_ultrasound_routine_anc_hm': 0,
                                                                                                       'preeclampsia_risk_screen___prevention_drug_routine_preecl_hm':preeclampsia_screening_prevention_coverage,
                                                                                                       'noninvasive_hb_diagnostic_first_trimester_an_hm':novel_first_trimester_hb_test_scale_up, 'noninvasive_hb_diagnostic_second_trimester_an_hm': novel_second_trimester_hb_test_scale_up,
                                                                                                       'ai_enabled_ultrasound_routine_anc_hm': ai_ultrasound_scale_up
                                                                                                       })

    return {
        'instructions_scaleup': instructions_scaleup,
        'instructions_scaleup_non_diagnosis': instructions_scaleup_non_diagnosis,
        'instructions_scaleup_diagnosis_only': instructions_scaleup_diagnosis_only,
    }
