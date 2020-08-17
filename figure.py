import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#from matplotlib import font_manager
manager = matplotlib.font_manager.FontManager()

# font_dirs = ['/extdata6/Minsu/bin/font']
# font_files = matplotlib.font_manager.findSystemFonts(fontpaths=font_dirs)
# for font_path in font_files:
#     manager.addfont(path=font_path)

#manager.addfont('/extdata6/Minsu/bin/font/Helvetica/Helvetica.ttf')
#manager.addfont('/extdata6/Minsu/bin/font/HMFMPYUN.ttf')
# path = '/extdata6/Minsu/bin/font/Helvetica/Helvetica.ttf'
# font = matplotlib.ft2font.FT2Font(path)
# helvetica_entry = matplotlib.font_manager.ttfFontProperty(font)
import sys
# print(helvetica_entry.name  )
# finding_property = matplotlib.font_manager.FontProperties(family='Helvetica')
# print(finding_property.get_family(), finding_property.get_file())
# print ( manager.score_family( finding_property.get_family(),  helvetica_entry.name   ) )
# print( manager.score_style(finding_property.get_style(), helvetica_entry.style))
# print( manager.score_variant(finding_property.get_variant(), helvetica_entry.variant))
# print( manager.score_weight(finding_property.get_weight(), helvetica_entry.weight))
# print( manager.score_stretch(finding_property.get_stretch(), helvetica_entry.stretch))
# print( manager.score_size(finding_property.get_size(), helvetica_entry.size))
# matplotlib.rc( 'text', usetex=False )
matplotlib.rc( 'font', family='Helvetica' )
# matplotlib.rc( 'font', family=manager.ttflist[-1].name )
# matplotlib.rc( 'font', size=30)
# matplotlib.rc( 'ps', useafm=True)
# matplotlib.rc( 'savefig', pad_inches = 0.2)

# manager.ttflist[-1].name
# print(manager.findfont(finding_property))# , fallback_to_default=False))
# default_prop =  matplotlib.font_manager.FontProperties()
# print(default_prop)

# print ( manager.score_family( default_prop.get_family(),  helvetica_entry.name   ) )
# print( manager.score_style(default_prop.get_style(), helvetica_entry.style))
# print( manager.score_variant(default_prop.get_variant(), helvetica_entry.variant))
# print( manager.score_weight(default_prop.get_weight(), helvetica_entry.weight))
# print( manager.score_stretch(default_prop.get_stretch(), helvetica_entry.stretch))
# print( manager.score_size(default_prop.get_size(), helvetica_entry.size))


#Helvetica_propety = font_manager.FontProperties(fname='/extdata6/Minsu/bin/font/Helvetica/Helvetica.ttf')
#helvetica_entry = matplotlib.font_manager.FontEntry( fname='/extdata6/Minsu/bin/font/Helvetica/Helvetica.ttf', )

#manager.ttflist.append( )
#print( Helvetica_propety.get_family() )

# print( font_manager.findfont(Helvetica_propety, directory='/extdata6/Minsu/bin/font/', fallback_to_default=False ) )
# FM.addfont(path='/extdata6/Minsu/bin/font/Helvetica/Helvetica.ttf')
# #font_manager.json_dump(FM, '/extdata6/Minsu/bin/font/font_list.json')
# font_manager.json_load('/extdata6/Minsu/bin/font/font_list.json')

# matplotlib.rc( 'font', family='sans-serif' )
#matplotlib.rc( 'font', {family='Helvetica' )


samples = ['AGS5', 'AGS1', 'AGS2', 'AGS4']
sample_dict = {'control':'Control', 'AGS5':'AGS5 (SAMHD1)', 'AGS1':'AGS1 (TREX1)', 'AGS2':'AGS2 (RNASEH2B)', 'AGS4':'AGS4 (RNASEH2A)'}
replicates = ('P1', 'P2')

read_counts = {
    'control' : 1.7465040,
    'AGS1'    : 2.1518048,
    'AGS2'    : 1.1997181,
    'AGS4'    : 2.2166405,
    'AGS5'    : 3.0926593,
}

color_dict = { 
    'control':'#F5B50B',
    'AGS1':'#C40234',
    'AGS2':'#87C756',
    'AGS4':'#8C4F92',
    'AGS5':'#3A73C0',
}

def load_drip_signal_array(sample_name) :
    file = './drip-table/%s_drip.tsv' % sample_name

    out_list = []
    for line in open(file) :
        signal_for_one_site =  [float(i) for i in line.split('\t') ]
        assert len(signal_for_one_site) == 240
        out_list.append( signal_for_one_site )
    return out_list

drip_data = {}
gene_center_data = {}
random_data = {}
for sample in samples+['control'] :
    drip_data[sample]        = np.concatenate([np.array(load_drip_signal_array(sample+'_CD_HO'                ))[:, ::-1]/read_counts[sample],
                                               np.array(load_drip_signal_array(sample+'_HO_CD'                ))         /read_counts[sample]])
    gene_center_data[sample] = np.concatenate([np.array(load_drip_signal_array(sample+'_gene_center_antisense'))[:, ::-1]/read_counts[sample],
                                               np.array(load_drip_signal_array(sample+'_gene_center_sense'    ))         /read_counts[sample]])
    random_data[sample]      = np.concatenate([np.array(load_drip_signal_array(sample+'_random_part1'         ))         /read_counts[sample], 
                                               np.array(load_drip_signal_array(sample+'_random_part2'         ))         /read_counts[sample], 
                                               np.array(load_drip_signal_array(sample+'_random_part3'         ))         /read_counts[sample]])
    ## random dataset are seperated to 3 files because of file-size limitation(25MB)

def bootstrap_one_sample(variable_1, n_repeats=10000, verbose=False):
    resampled_var1 = np.zeros((n_repeats, variable_1.shape[1]), dtype=np.float32)
    for i in range(n_repeats):
        if verbose and not i % (n_repeats//10): ## progress check\
            print(i)
        random_sample = np.random.choice(len(variable_1) - 1, len(variable_1) - 1, replace=True)
        resampled_var1[i, :] = variable_1[random_sample, :].mean(axis=0)
    return resampled_var1

def plot_mean_with_95ci(bootstrap_sample, color, ax, x=None):
    nrow, ncol = bootstrap_sample.shape    
    if x is None:
        x = np.linspace(-12, 12, ncol)
    # lower_bound = np.percentile(bootstrap_sample, 2.5, axis=0)
    # upper_bound = np.percentile(bootstrap_sample, 97.5, axis=0)
    # ax.fill_between(x, lower_bound, upper_bound, color=color, linestyle='None', alpha=0.3)
    ax.plot(x, bootstrap_sample.mean(axis=0), color=color)

drip_boot = { sample : bootstrap_one_sample(drip_data[sample]) for sample in samples+['control'] }
gene_center_boot = { sample : bootstrap_one_sample(gene_center_data[sample]) for sample in samples+['control'] }
random_boot = { sample : bootstrap_one_sample(random_data[sample]) for sample in samples+['control'] }

def figure1a(drip_boot, gene_center_boot, random_boot):
    image_count = 0
    fig, ax = plt.subplots(nrows=4, ncols=3, sharex=True, figsize=(22,20))

    y_axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=None, useMathText=True, useLocale=None)
    y_axis_formatter.set_powerlimits( (-1,1) )
    y_axis_formatter.set_scientific( True )
    ylim_range = (0.24, 0.55)

    for sample in samples :
        ## 1st column : Replication origin
        plot_mean_with_95ci(drip_boot['control'], color_dict['control'], ax[image_count][0])
        plot_mean_with_95ci(drip_boot[sample]   , color_dict[sample], ax[image_count][0])

        ax[image_count][0].axvline(x=0, linestyle='dotted', color='k', linewidth=0.5)

        ax[image_count][0].set_ylabel('average DRIP-seq readcount', fontsize=12)
        ax[image_count][0].yaxis.set_major_formatter(y_axis_formatter)
        ax[image_count][0].yaxis.set_offset_position( 'left' )
        ax[image_count][0].yaxis.tick_right()
        ax[image_count][0].set_ylim(ylim_range)
        ax[image_count][0].set_xlim(-12,12)
        
        if image_count == 3 : ax[image_count][0].set_xlabel('Distance from replication origin (kb)', fontsize=12)


        ## 2nd column : Center of Gene bodies
        plot_mean_with_95ci(gene_center_boot['control'], color_dict['control'], ax[image_count][1])
        plot_mean_with_95ci(gene_center_boot[sample]   , color_dict[sample], ax[image_count][1])

        ax[image_count][1].axvline(x=0 , linestyle='dotted', color='k', linewidth=0.5)
        
        ax[image_count][1].yaxis.set_major_formatter(y_axis_formatter)
        ax[image_count][1].yaxis.set_offset_position( 'left' )
        ax[image_count][1].yaxis.tick_right()
        ax[image_count][1].set_ylim(ylim_range)
        
        ax[image_count][1].set_xlim(-12,12)
        if image_count == 3 : ax[image_count][1].set_xlabel('Distance from center of gene bodies (kb)', fontsize=12)


        ## 3rd column : Center of Gene bodies
        plot_mean_with_95ci(random_boot['control'], color_dict['control'], ax[image_count][2])
        plot_mean_with_95ci(random_boot[sample]   , color_dict[sample], ax[image_count][2])

        ax[image_count][2].axvline(x=0 , linestyle='dotted', color='k', linewidth=0.5)
        
        ax[image_count][2].yaxis.set_major_formatter(y_axis_formatter)
        ax[image_count][2].yaxis.set_offset_position( 'left' )
        ax[image_count][2].yaxis.tick_right()
        ax[image_count][2].set_ylim(ylim_range)
        
        ax[image_count][2].set_xlim(-12,12)
        if image_count == 3 : ax[image_count][2].set_xlabel('Distance from random position (kb)', fontsize=12)


        ## legend
        figure_handles = [matplotlib.patches.Patch(color=color_dict['control'], alpha=0.3),
                          matplotlib.patches.Patch(color=color_dict[sample], alpha=0.3)]
        ax[image_count][2].legend(labels=['Control', sample_dict[sample] ], bbox_to_anchor=(1.1, 0.99), handles = figure_handles, loc='upper left', fontsize=12)

        image_count +=1 

    plt.tight_layout()
    
    out_fig_route = os.path.join('.','figure','fig1a.png')
    plt.savefig(out_fig_route)
    
    plt.close()

def figure1b(drip_boot) : 
    import matplotlib.cbook as cbook
    graph_data = []
    graph_ho_data =[]
    graph_cd_data = []
    ylim_range = (0.29, 0.51)
    for sample in ['control']+samples : 
        nrow, ncol = drip_boot[sample].shape
        sample_data = []
        sample_ho_data = []
        sample_cd_data = []
        for i in range(nrow) :
            assert len(drip_boot[sample][i, :]) == ncol, len(drip_boot[sample][i, :])
            sample_data.append( (drip_boot[sample][i, 60:180]).mean())
            sample_ho_data.append( (drip_boot[sample][i, 60:120]).mean())
            sample_cd_data.append( (drip_boot[sample][i, 120:180]).mean())
        assert len(sample_data) == nrow
        stat_data = cbook.boxplot_stats(sample_data)[0]
        graph_data.append(stat_data)

        stat_ho_data = cbook.boxplot_stats(sample_ho_data)[0]
        graph_ho_data.append(stat_ho_data)

        stat_cd_data = cbook.boxplot_stats(sample_cd_data)[0]
        graph_cd_data.append(stat_cd_data)

    y_axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=True, useMathText=True, useLocale=None)
    y_axis_formatter.set_powerlimits( (-1,1) )
    y_axis_formatter.set_scientific( True )

    ax = plt.axes()
    ax.bxp(graph_data, widths=0.3, showfliers=False)
    plt.xticks([1,2,3,4,5], [ sample_dict[i].replace(' ', '\n') for i in ['control']+samples], fontsize=12)
    plt.title('Average DRIP-seq readcount in 12kb window', fontsize=12)
    plt.ylabel('average DRIP-seq readcount', fontsize=12)
    plt.ylim(ylim_range)
    ax.yaxis.set_major_formatter(y_axis_formatter)
    ax.yaxis.set_offset_position( 'left' )
    out_fig_route = os.path.join('.','figure','fig1b.png')
    plt.savefig(out_fig_route)
    plt.close()

    ax = plt.axes()
    ax.bxp(graph_ho_data, widths=0.3, showfliers=False)
    plt.xticks([1,2,3,4,5], [ sample_dict[i].replace(' ', '\n') for i in ['control']+samples], fontsize=12)
    plt.title('Average DRIP-seq readcount in 6kb HO region', fontsize=12)
    plt.ylabel('average DRIP-seq readcount', fontsize=12)
    plt.ylim(ylim_range)
    ax.yaxis.set_major_formatter(y_axis_formatter)
    ax.yaxis.set_offset_position( 'left' )
    out_fig_route = os.path.join('.','figure','fig1b_ho.png')
    plt.savefig(out_fig_route)
    plt.close()

    ax = plt.axes()
    ax.bxp(graph_cd_data, widths=0.3, showfliers=False)
    plt.xticks([1,2,3,4,5], [ sample_dict[i].replace(' ', '\n') for i in ['control']+samples], fontsize=12)
    plt.title('Average DRIP-seq readcount in 6kb CD region', fontsize=12)
    plt.ylabel('average DRIP-seq readcount', fontsize=12)
    plt.ylim(ylim_range)
    ax.yaxis.set_major_formatter(y_axis_formatter)
    ax.yaxis.set_offset_position( 'left' )
    out_fig_route = os.path.join('.','figure','fig1b_cd.png')
    plt.savefig(out_fig_route)
    plt.close()

def figure1c(drip_data):
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    for sample in ['control']+samples : 
        nrow, ncol = drip_data[sample].shape
        sample_data = []
        for i in range(nrow) :
            #assert len(drip_data[sample][: , i]) == nrow, len(drip_data[sample][i, :])
            sample_data.append( (drip_data[sample][i, 60:180] ).sum()  )

        graph_data = []
        bin_set = list( set(sample_data) )
        bin_set.sort()
        for bins in bin_set :
            count = 0 
            for data in sample_data :
                if data < bins : 
                    count+=1
            graph_data.append(count/len(sample_data))
        ax.step( bin_set, graph_data, color=color_dict[sample], label=sample)

    ax.legend(loc='lower right')

    ax.set_ylabel('cumulative distribution frequency')
    ax.set_xlabel('mapped read in 12kb windows for each origins')

    ax.set_xlim(0,200)

    ax.set_ylim(0,1)
    # x_axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=0.01, useLocale=None)
    # x_axis_formatter.set_powerlimits( (-2, -2) )
    # x_axis_formatter.set_scientific( True )
    # ax.xaxis.set_major_formatter(x_axis_formatter)

    fig.tight_layout()
    out_fig_route = os.path.join('.','figure','fig1c.png')
    fig.savefig(out_fig_route)

figure1a(drip_boot, gene_center_boot, random_boot)
figure1b(drip_boot)
figure1c(drip_data)
