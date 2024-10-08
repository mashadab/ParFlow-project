from parflow.tools.io import write_pfb, read_pfb

satur_data = read_pfb("sat_fname")


press_data = sat2press(satur_data)



write_pfb(press_data, "ic_pressure.pfb")