import atomica as at

def fill_assumptions_as_data(P, years=[2020]):
    at.logger.info(
        f'Adding "data" values for years {years} that match the assumption for that parameter for all parameters with an assumption and no data in project {P.name}')
    # replace in parsets
    for parset in P.parsets.keys():
        for par in P.parsets[parset].pars.keys():
            for pop in P.parsets[parset].pars[par].ts.keys():
                if P.parsets[parset].pars[par].ts[pop].has_data and not P.parsets[parset].pars[par].ts[
                    pop].has_time_data:
                    P.parsets[parset].pars[par].ts[pop].t = list(years)
                    P.parsets[parset].pars[par].ts[pop].vals = [P.parsets[parset].pars[par].ts[pop].assumption for _ in
                                                                years]

    # replace in data as well
    for par in P.data.tdve.keys():
        for pop in P.data.tdve[par].ts.keys():
            if P.data.tdve[par].ts[pop].has_data and not P.data.tdve[par].ts[pop].has_time_data:
                P.data.tdve[par].ts[pop].t = list(years)
                P.data.tdve[par].ts[pop].vals = [P.data.tdve[par].ts[pop].assumption for _ in years]

    return P