file = open('obs_grouped.txt', 'r')
group = []
mos = ['jan', 'feb', 'mar', 'apr', 'may', 'jul', 'aug', 'sep', 'oct', 'nov', 'jun', 'dec']

c = 0
for line in file:
    c += 1
    print(c)
    if line == '\n':
        for mo in mos:
            if mo in group[0]:
                month = mo
                year = group[0][group[0].find(mo)-2:group[0].find(mo)].strip('/')
                date = group[0][group[0].find(mo):group[0].find(mo)+6].replace('/', '').replace(mo, '')
                for letter in 'qwertyuiopazsdfghjklzxcvbnm':
                    date = date.replace(letter, '')
                break
        wavelens = {}
        for obs in group:
            wavelen = float(obs[obs.find('2.'):obs.find('2.')+4])
            if wavelen in wavelens:
                wavelens[wavelen].append(obs)
            else:
                wavelens[wavelen] = [obs]
        wavelens['misc'] = []
        delwavelens = []
        for wavelen in wavelens.keys():
            if len(wavelens[wavelen]) < 3:
                for obs in wavelens[wavelen]:
                    wavelens['misc'].append(obs)
                delwavelens.append(wavelen)
        for wavelen in delwavelens:
            del wavelens[wavelen]
        for wavelen in wavelens.keys():
            filename = f'{year}-{month}-{date}_{wavelen}.txt'
            outfile = open(filename, 'w')
            for obs in wavelens[wavelen]:
                outfile.write(obs+'\n')
            outfile.close()
        group = []
    else:
        group.append(line.strip('\n'))
