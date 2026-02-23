## PLACEHOLDER - this script should create all of the databooks

if not databook.exists():
    pops = {
        'anc_hm': {'label': 'Antenatal care: Humans', 'type': 'anc_hm'},
        'ol_hm': {'label': 'Obstructed labor: Humans', 'type': 'ol_hm'},
        'pt_hm': {'label': 'Preterm/SGA: Humans', 'type': 'pt_hm'},
        'msep_hm': {'label': 'Maternal sepsis: Humans', 'type': 'msep_hm'},
        'hm_hm': {'label': 'Hemorrhage: Humans', 'type': 'hm_hm'},
        'asx_hm': {'label': 'Asphyxia: Humans', 'type': 'asx_hm'},
        'oth_hm': {'label': 'Other neonatal conditions: Humans', 'type': 'oth_hm'},
        'nsep_hm': {'label': 'Neonatal sepsis: Humans', 'type': 'nsep_hm'},
        'an_hm': {'label': 'Anemia: Humans', 'type': 'an_hm'},
        'muw_hm': {'label': 'Maternal underweight: Humans', 'type': 'muw_hm'},
        'wasting_hm': {'label': 'Wasting: Humans', 'type': 'wasting_hm'},
        'stunting_hm': {'label': 'Stunting: Humans', 'type': 'stunting_hm'},
        'rds_hm': {'label': 'Preterm/SGA: Humans', 'type': 'rds_hm'},
    }
    D = at.ProjectData.new(F, np.arange(2020, 2040), pops=pops, transfers=0)
    D.save(databook)
