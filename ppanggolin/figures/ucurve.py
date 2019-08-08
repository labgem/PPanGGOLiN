#default libraries
import logging
from collections import defaultdict

#installed libraries
import plotly.graph_objs as go
import plotly.offline as out_plotly
#local libraries
from ppanggolin.formats import checkPangenomeInfo

def drawUCurve(pangenome, output, soft_core = 0.95):
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=True)
    logging.getLogger().info("Drawing the U-shaped curve...")
    max_bar = 0
    count = defaultdict(lambda : defaultdict(int))
    is_partitionned = False
    has_undefined = False
    for fam in pangenome.geneFamilies:
        nb_org  = len(fam.organisms)
        if fam.partition != "":
            is_partitionned = True
            if fam.partition == "U":
                has_undefined = True
            count[nb_org][fam.namedPartition]+=1
        count[nb_org]["pangenome"]+=1
        max_bar = count[nb_org]["pangenome"] if count[nb_org]["pangenome"] > max_bar else max_bar
    data_plot = []
    chao = "NA"
    if count[1]["pangenome"] > 0:
        chao = round(len(pangenome.geneFamilies) + ((count[0]["pangenome"]^2)/(count[1]["pangenome"]*2)),2)
    COLORS = {"pangenome":"black", "exact_accessory":"#EB37ED", "exact_core" :"#FF2828", "soft_core":"#c7c938", "soft_accessory":"#996633","shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}
    if is_partitionned and not has_undefined:
        persistent_values = []
        shell_values      = []
        cloud_values      = []
        for nb_org in range(1,len(pangenome.organisms)+1):
            persistent_values.append(count[nb_org]["persistent"])
            shell_values.append(count[nb_org]["shell"])
            cloud_values.append(count[nb_org]["cloud"])
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=persistent_values,name='persistent', marker=dict(color = COLORS["persistent"])))
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=shell_values,name='shell', marker=dict(color = COLORS["shell"])))
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=cloud_values,name='cloud', marker=dict(color = COLORS["cloud"])))
    else:
        text = 'undefined' if has_undefined else "pangenome"
        undefined_values = []
        for nb_org in range(1,len(pangenome.organisms)+1):
            undefined_values.append(count[nb_org][text])
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=undefined_values,name=text, marker=dict(color = COLORS[text])))
    layout = None
    x = len(pangenome.organisms)*soft_core
    layout =  go.Layout(title = "Gene families frequency distribution (U shape), chao="+str(chao),
                        xaxis = dict(title='Occurring in x genomes'),
                        yaxis = dict(title='# of gene families (F)'),
                        barmode='stack', shapes=[dict(type='line', x0=x, x1=x, y0=0, y1=max_bar, line = dict(dict(width=5, dash='dashdot', color="grey")))])

    fig = go.Figure(data=data_plot, layout=layout)
    out_plotly.plot(fig, filename = output+"/Ushaped_plot.html", auto_open=False)
    logging.getLogger().info(f"Done drawing the U-shaped curve : '{output+'/Ushaped_plot.html'}'")
