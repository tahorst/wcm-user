
subunits50s = sim_data.process.complexation.getMonomers("CPLX0-3962[c]")["subunitIds"]
idx = [np.where(sim_data.process.translation.monomerData["id"] == x)[0][0] for x in subunits50s if x in sim_data.process.translation.monomerData["id"]]
aa50s = sim_data.process.translation.monomerData[idx]["aaCounts"].asNumber()
aa50s[-2, :] *= 2

subunits30s = sim_data.process.complexation.getMonomers("CPLX0-3953[c]")["subunitIds"]
idx = [np.where(sim_data.process.translation.monomerData["id"] == x)[0][0] for x in subunits30s if x in sim_data.process.translation.monomerData["id"]]
aa30s = sim_data.process.translation.monomerData[idx]["aaCounts"].asNumber()

aaUsed = np.sum(aa50s, axis = 0) + np.sum(aa30s, axis = 0)
aaFraction = aaUsed / float(np.sum(aaUsed))
