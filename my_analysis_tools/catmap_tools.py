
def create_catmap_input(home,element,facet,scaling_data,basename='/catmap/catmap_input'):
    out= open(home+basename+'_%s.txt'%facet,'w')
    out.writelines('surface_name\tsite_name\tspecies_name\tformation_energy\tfrequencies\treference\n')
    out.writelines('None\tgas\tCO\t0.0\t[2136.4]\town calc\n')
    out.writelines('None\tgas\tpe\t0.0\t[]\town calc\n')
    #out.writelines('None\tgas\tCO\t2.5316863115494925\t[2152.5]\town calc\n')
    out.writelines('None\tgas\tH2\t0.0\t[4426.5]\town calc\n')
    out.writelines('None\tgas\tH2O\t0.0\t[1616.2,3741.5,3855.6]\town calc\n')
    #out.writelines('None\tgas\tCH4\t0.0\t[1327.7, 1329.8, 1332.2, 1546.0, 1547.6, 2983.5, 3118.2, 3120.6, 3121.4]\town calc\n')
    out.writelines('None\tgas\tC2H4\t-2.76975805\t[826.36096307,  959.45228605,  960.34934469, 1053.97242644, '
             '1229.69837936, 1348.85228724, 1472.24547713, 1668.96347041, 3105.59055667, '
             '3125.65012816, 3173.82607546, 3196.77938303]\town calc\n')
    out.writelines('None\tgas\tCH4\t-2.5162865\t[]\town calc\n')

    for element in scaling_data.keys():
        #out.writelines('%s\t%s\tCO\t%s\t[]\town calc\n'
        for ads in scaling_data[element].keys():
          if ads not in ['CO','OCCO','OCHCO','OCCOH','H','CCO','HCCO','OC-CO','OCCHOH','CHO','COH']:
            out.writelines('%s\t%s\t%s\t%s\t[]\town calc\n'
                       %(element,facet,ads,scaling_data[element][ads]))

        if 'CO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tCO\t%s\t[15.5,18.8,251.3,252.3,303.1,1952.8]\town calc\n'
                       %(element,facet,scaling_data[element]['CO']))
        if 'OCCO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tOCCO\t%s\t[86.7,109.3,115.1,218.0,227.2,243.9,317.6,321.0,552.7,645.6,1407.3,1451.8]\town calc\n'
                       %(element,facet,scaling_data[element]['OCCO']))
        if 'OCHCO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tOCHCO\t%s\t[7.0,29.4,42.4,69.5,188.4,224.5,249.0,410.5,638.0,806.2,914.2,1312.3,1620.8,1714.8,2906.4]\town calc\n'
                       %(element,facet,scaling_data[element]['OCHCO']))
        if 'OCCOH' in scaling_data[element].keys():
            out.writelines('%s\t%s\tOCCOH\t%s\t[93.4,108.5,145.2,218.6,256.0,271.9,322.0,373.5,493.2,658.2,811.6,1021.0,1252.2,1543.1,3691.1]\town calc\n'
                       %(element,facet,scaling_data[element]['OCCOH']))
        if 'H' in scaling_data[element].keys():
            out.writelines('%s\t%s\tH\t%s\t[307.7,308.8,636.3]\town calc\n'
                       %(element,facet,scaling_data[element]['H']))
        if 'CCO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tCCO\t%s\t[49.0,51.6,202.8,203.5,496.1,496.2,671.9,1172.0,1908.6]\town calc\n'
                       %(element,facet,scaling_data[element]['CCO']))
        if 'HCCO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tHCCO\t%s\t[48.0, 49.7,109.4,231.3,305.7,508.3,551.4,565.6,970.7,1165.3,2008.0,3099.1]\town calc\n'
                       %(element,facet,scaling_data[element]['HCCO']))
        if 'OC-CO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tOC-CO\t%s\t[90.8,101.4,115.6,214.2,222.1,230.4,250.1,302.0,516.2,1568.8,1599.3]\town calc\n'
                       %(element,facet,scaling_data[element]['OC-CO']))
        if 'OCCHOH' in scaling_data[element].keys():
            out.writelines('%s\t%s\tOCCHOH\t%s\t[12.0,75.8,182.0,211.4,222.6,230.3,328.0,468.8,547.2,756.9,811.4,918.1,1048.5,1231.9,1303.1,1456.0,3161.2,3650.8]\town calc\n'
                       %(element,facet,scaling_data[element]['OCCHOH']))
        if 'CHO' in scaling_data[element].keys():
            out.writelines('%s\t%s\tCHO\t%s\t[18.0,32.8,76.9,190.3,400.4,646.7,1227.0,1581.5,2729.7]\town calc\n'
                       %(element,facet,scaling_data[element]['CHO']))
        if 'COH' in scaling_data[element].keys():
            out.writelines('%s\t%s\tCOH\t%s\t[15,139.8,159.6,273.0,366.9,402.2,951.4,1171.8,3648.3]\town calc\n'
                       %(element,facet,scaling_data[element]['COH']))

        #Careful OCCHOH and COH had an imaginary


    out.close()
