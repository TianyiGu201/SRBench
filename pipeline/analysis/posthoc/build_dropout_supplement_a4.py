"""Render dropout curves on A4 portrait pages.

The complete 3-row × 2-column method group is moved vertically as one unit
using its measured artist bounding box. No significance labels appear in figures.
Run after build_dropout_paired_statistics.py. The default raw mode corresponds to the final manuscript figures; delta mode is an optional descriptive view.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED_CSV = HERE / 'dropout_seed_paired_effects.csv'
TABLE_FILE = HERE / 'Gu_Supplementary_Table_6.xlsx'
AXIS_FILE = HERE / 'dropout_final_axis_limits.csv'
COMBINED_FILE = HERE / 'Gu_Supplementary_Figure_Dropout_A4.pdf'
MAIN_FILE = HERE / 'Gu_Main_Figure_5d-i_Dropout_Revision.pdf'
SUPP_FILE = HERE / 'Gu_Supplementary_Figures_Dropout_Datasets2-6.pdf'
CROP_FILE = HERE / 'Gu_Main_Figure_5d-i_Dropout_Revision_Cropped.pdf'
AUDIT_FILE = HERE / 'dropout_figure_layout_validation.csv'
DATASETS = [('Xe_hsa_breast_c',1),('Xe_hsa_colon_h',2),('Xe_hsa_kidney_h',3),
            ('Xe_hsa_heart_h',4),('Xe_hsa_skin_h',5),('Xe_hsa_liver_h',6)]
METHODS = ['XFuse','TESLA','iSTAR','scstGCN','SuperST_t','SuperST_t_g']
METRICS = ['PCC','SCC','JSD','KLD']
RATIOS = [0.1,0.2,0.3,0.4,0.5]
COLORS = {'PCC':'#F8766D','SCC':'#20B2AA','JSD':'#7570B3','KLD':'#E69F00'}
MM_TO_IN = 1/25.4
PANEL_FONT = FontProperties(family='Arial',weight='bold',size=8)
plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':6,
    'axes.labelsize':6,'axes.titlesize':6,'xtick.labelsize':6,'ytick.labelsize':6,
    'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.55})


def parse_ticks(value):
    labels = [part.strip() for part in value.split(',')]
    return [float(part) for part in labels], labels


def raw_metric_scales(raw):
    """Derive one readable raw-metric scale for every dataset × metric."""
    records=[]
    for (dataset,metric), group in raw.groupby(['Dataset_ID','Metric'], sort=False):
        lo=float(group.Value.min())
        hi=float(group.Value.max())
        spread=hi-lo
        pad=0.05*spread if spread>0 else max(abs(lo)*0.05,0.01)
        locator=MaxNLocator(nbins=4,steps=[1,2,2.5,5,10],min_n_ticks=3)
        ticks=locator.tick_values(lo-pad,hi+pad)
        if len(ticks)<2:
            ticks=np.linspace(lo-pad,hi+pad,4)
        ymin=float(ticks[0])
        ymax=float(ticks[-1])
        # Two decimals retain the readable precision needed for these compact
        # panels while preventing tick labels from encroaching on neighbors.
        decimals=2 if np.max(np.abs(ticks))<1 else 1
        labels=[f'{tick:.{decimals}f}' for tick in ticks]
        records.append({'dataset':dataset,'metric':metric,'suggested_ymin':ymin,
                        'suggested_ymax':ymax,'visible_ticks':','.join(map(str,ticks)),
                        'tick_labels':','.join(labels)})
    return pd.DataFrame(records)


def artist_bbox(artist, renderer):
    return artist.get_tightbbox(renderer) if hasattr(artist,'get_tightbbox') else artist.get_window_extent(renderer)


def block_bbox(fig, renderer, rows=None):
    if rows is None:
        artists = fig.axes + fig.texts
    else:
        artists = [ax for ax in fig.axes if ax._method_row in rows]
        artists += [tx for tx in fig.texts if tx._method_row in rows]
    boxes = [artist_bbox(a,renderer) for a in artists if a.get_visible()]
    return Bbox.union(boxes)


def shift_composite(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    before = block_bbox(fig,renderer)
    dy_normalized = (fig.bbox.height/2 - before.y0/2 - before.y1/2)/fig.bbox.height
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0,pos.y0+dy_normalized,pos.width,pos.height])
    for tx in fig.texts:
        x,y = tx.get_position()
        tx.set_position((x,y+dy_normalized))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    after = block_bbox(fig,renderer)
    page_center_error_pt = (after.y0+after.y1-fig.bbox.height)/2 * 72/fig.dpi
    top_margin_pt = (fig.bbox.height-after.y1)*72/fig.dpi
    bottom_margin_pt = after.y0*72/fig.dpi
    row_boxes = [block_bbox(fig,renderer,{row}) for row in range(3)]
    row_gaps_pt = [(row_boxes[i].y0-row_boxes[i+1].y1)*72/fig.dpi for i in (0,1)]
    assert abs(page_center_error_pt) < 1.0, page_center_error_pt
    assert abs(top_margin_pt-bottom_margin_pt) < 2.0
    assert min(row_gaps_pt) > 12.0, row_gaps_pt
    assert after.x0 > 0 and after.x1 < fig.bbox.width and after.y0 > 0 and after.y1 < fig.bbox.height
    return after, {'center_error_pt':page_center_error_pt,'top_margin_pt':top_margin_pt,
                    'bottom_margin_pt':bottom_margin_pt,'row1_row2_gap_pt':row_gaps_pt[0],
                    'row2_row3_gap_pt':row_gaps_pt[1],
                    'block_height_pt':after.height*72/fig.dpi}


def render_page(dataset, number, seeds, summary, scales, raw_metrics=False):
    fig = plt.figure(figsize=(210*MM_TO_IN,297*MM_TO_IN),facecolor='white')
    width_pt,height_pt = fig.get_size_inches()*72
    # A4 portrait artwork must fit between Illustrator guides x=18 and 192 mm.
    # Retain the 6/7/8 pt fonts and square metric areas; compress geometry only.
    group_lefts_pt=(50.6,296.6)
    group_tops_pt=(785.0,520.0,255.0)
    plot_side_pt=33.0
    plot_pitch_pt=62.5
    panel_letters='defghi' if number==1 else 'abcdef'
    audit=[]
    for method_idx, method in enumerate(METHODS):
        row,col=divmod(method_idx,2)
        top_pt=group_tops_pt[row]
        left_group=group_lefts_pt[col]
        letter=fig.text((left_group+2)/width_pt,top_pt/height_pt,panel_letters[method_idx],
            ha='left',va='top',fontproperties=PANEL_FONT,color='black')
        letter._method_row=row
        title=fig.text((left_group+135.5)/width_pt,top_pt/height_pt,
            f'Dropout-Associated Prediction Performance of\n{method} in Dataset {number}',
            ha='center',va='top',fontsize=7,fontweight='normal',linespacing=1.02,color='black')
        title._method_row=row
        for metric_idx,metric in enumerate(METRICS):
            left_pt=left_group+25+metric_idx*plot_pitch_pt
            scale_row=scales[(scales.dataset==dataset)&(scales.metric==metric)]
            assert len(scale_row)==1
            scale=scale_row.iloc[0]
            yticks,_=parse_ticks(scale.visible_ticks)
            yticklabels=scale.tick_labels.split(',') if 'tick_labels' in scale else parse_ticks(scale.visible_ticks)[1]
            sub=summary[(summary.Dataset_ID==dataset)&(summary.Method==method)&(summary.Metric==metric)]
            ax=fig.add_axes([left_pt/width_pt,(top_pt-26-plot_side_pt)/height_pt,
                             plot_side_pt/width_pt,plot_side_pt/height_pt])
            ax._method_row=row
            ax.set_xlim(0.075,0.525)
            ax.set_ylim(float(scale.suggested_ymin),float(scale.suggested_ymax))
            ax.set_xticks(RATIOS,[f'{r:.1f}' for r in RATIOS])
            ax.set_yticks(yticks,yticklabels)
            ax.set_xlabel('Dropout ratio',fontsize=6,labelpad=1.0,fontweight='normal',color='black')
            ax.set_ylabel(metric if raw_metrics else f'Δ{metric}',fontsize=6,labelpad=0.3,fontweight='normal',color='black')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(axis='both',direction='out',length=1.2,width=0.5,pad=1.0,
                           colors='black',labelsize=6)
            ax.tick_params(axis='x',labelrotation=45)
            for label in ax.get_xticklabels():
                label.set_horizontalalignment('right')
                label.set_fontweight('normal')
            for label in ax.get_yticklabels():
                label.set_fontweight('normal')
            ax.grid(False)
            if not raw_metrics:
                ax.axhline(0,color='#666666',linewidth=0.55,zorder=1)
            if len(sub)==5 and (raw_metrics or (sub.Status=='Available').all()):
                assert len(sub)==5
                sub=sub.set_index('Dropout_ratio').loc[RATIOS]
                for seed in (123,124,125):
                    vals=seeds[(seeds.Dataset_ID==dataset)&(seeds.Method==method)
                        &(seeds.Metric==metric)&(seeds.Seed==seed)].set_index('Dropout_ratio').loc[RATIOS,'Value'].to_numpy()
                    assert np.isfinite(vals).all()
                    ax.plot(RATIOS,vals,'o-',color='#929292',alpha=0.48,
                        markersize=1.5,linewidth=0.50,zorder=2)
                means=sub.Mean.to_numpy()
                sd=sub.SD.to_numpy()
                artist=ax.errorbar(RATIOS,means,yerr=sd,fmt='o-',color=COLORS[metric],
                    ecolor=COLORS[metric],linewidth=0.85,elinewidth=0.70,capsize=1.0,
                    capthick=0.55,markersize=2.1,zorder=4)
                assert np.array_equal(np.asarray(artist.lines[0].get_ydata(),dtype=float),means)
                observed=seeds[(seeds.Dataset_ID==dataset)&(seeds.Metric==metric)
                    &(seeds.Method==method)].Value
                assert observed.min()>=ax.get_ylim()[0] and observed.max()<=ax.get_ylim()[1]
                status='Available'
            else:
                assert number==5 and method=='scstGCN' and len(sub) in (0,5)
                ax.text(0.5,0.5,'Not\nevaluated',transform=ax.transAxes,ha='center',va='center',
                    fontsize=6,fontweight='normal',color='black')
                status='Not evaluated'
            audit.append({'Dataset':f'Dataset {number}','Method':method,'Metric':metric,
                'Ymin':ax.get_ylim()[0],'Ymax':ax.get_ylim()[1],
                'Plotted_mean_equal_table':True,'Status':status})
    # Shared limits across methods within each dataset and metric.
    audit_df=pd.DataFrame(audit)
    for metric,grp in audit_df.groupby('Metric'):
        assert len(grp)==6 and grp[['Ymin','Ymax']].drop_duplicates().shape[0]==1,(dataset,metric)
    bbox,layout=shift_composite(fig)
    left_mm=bbox.x0/fig.dpi*25.4
    right_mm=bbox.x1/fig.dpi*25.4
    assert left_mm>=18.0 and right_mm<=192.0, (number,left_mm,right_mm)
    layout['left_content_mm']=left_mm
    layout['right_content_mm']=right_mm
    layout['content_width_mm']=right_mm-left_mm
    # Validate the full tight bounds of neighboring metric panels, including
    # ylabel and tick-label extents, after fitting into the 174-mm guide span.
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    horizontal_gaps=[]
    for row in range(3):
        row_axes=fig.axes[row*8:(row+1)*8]
        for left_ax,right_ax in zip(row_axes,row_axes[1:]):
            horizontal_gaps.append((right_ax.get_tightbbox(renderer).x0-
                                    left_ax.get_tightbbox(renderer).x1)*72/fig.dpi)
    assert min(horizontal_gaps)>2.0, (number,min(horizontal_gaps))
    layout['min_subplot_horizontal_gap_pt']=min(horizontal_gaps)
    # No test text is allowed among any displayed figure labels.
    labels=[t.get_text() for t in fig.texts]
    for ax in fig.axes:
        labels.extend([t.get_text() for t in ax.texts])
        labels.extend([ax.get_xlabel(),ax.get_ylabel()])
    assert not any(re.search(r'Holm|p\s*[=<]|P\s*[=<]|permutation|significan',x,re.I) for x in labels)
    audit_df['Page_center_error_pt']=layout['center_error_pt']
    audit_df['Top_margin_pt']=layout['top_margin_pt']
    audit_df['Bottom_margin_pt']=layout['bottom_margin_pt']
    audit_df['Row1_to_row2_gap_pt']=layout['row1_row2_gap_pt']
    audit_df['Row2_to_row3_gap_pt']=layout['row2_row3_gap_pt']
    audit_df['Block_height_pt']=layout['block_height_pt']
    audit_df['Left_content_mm']=layout['left_content_mm']
    audit_df['Right_content_mm']=layout['right_content_mm']
    audit_df['Content_width_mm']=layout['content_width_mm']
    audit_df['Min_subplot_horizontal_gap_pt']=layout['min_subplot_horizontal_gap_pt']
    return fig,bbox,audit_df


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--preview-dir',type=Path,default=None)
    parser.add_argument('--mode',choices=['raw','delta'],default='raw')
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--axis-file',type=Path,help='Required only for delta mode')
    args=parser.parse_args()
    global HERE
    HERE=args.input_dir
    args.raw_metrics=args.mode=='raw'
    args.output_dir.mkdir(parents=True,exist_ok=True)
    output_dir=args.output_dir
    if not args.raw_metrics and args.axis_file is None:
        parser.error('--axis-file is required in delta mode')
    if args.raw_metrics:
        raw=pd.read_csv(HERE/'dropout_replicate_means_all_Xenium.csv',float_precision='round_trip')
        seeds=raw.rename(columns={'dataset':'Dataset_ID','method':'Method','metric':'Metric',
                                  'dropout_ratio':'Dropout_ratio','replicate_gene_mean':'Value',
                                  'replicate':'Replicate'}).copy()
        seeds['Seed']=seeds.Replicate.map({1:123,2:124,3:125})
        summary=(seeds.groupby(['Dataset_ID','Method','Metric','Dropout_ratio'],as_index=False)
                 .Value.agg(Mean='mean',SD=lambda x:x.std(ddof=1)))
        scales=raw_metric_scales(seeds)
        combined_file=output_dir/'Gu_Supplementary_Figure_Dropout_OriginalMetric_SharedYAxis_A4.pdf'
        main_file=output_dir/'Gu_Main_Figure_5d-i_Dropout_OriginalMetric_SharedYAxis.pdf'
        supp_file=output_dir/'Gu_Supplementary_Figures_Dropout_OriginalMetric_SharedYAxis_Datasets2-6.pdf'
        crop_file=output_dir/'Gu_Main_Figure_5d-i_Dropout_OriginalMetric_SharedYAxis_Cropped.pdf'
        audit_file=output_dir/'dropout_original_metric_shared_axis_layout_validation.csv'
        scales.to_csv(output_dir/'dropout_original_metric_shared_axis_limits.csv',index=False)
    else:
        seeds=pd.read_csv(HERE/SEED_CSV.name, float_precision='round_trip').rename(columns={'Delta':'Value'})
        reported=pd.read_excel(HERE/TABLE_FILE.name, sheet_name='All_Dropout_Effects')
        summary=reported.rename(columns={'Mean_delta':'Mean','SD_delta':'SD'})
        scales=pd.read_csv(args.axis_file)
        combined_file,main_file,supp_file,crop_file,audit_file=tuple(output_dir/p.name for p in (COMBINED_FILE,MAIN_FILE,SUPP_FILE,CROP_FILE,AUDIT_FILE))
    seeds['Dropout_ratio']=seeds['Dropout_ratio'].round(1)
    summary['Dropout_ratio']=summary['Dropout_ratio'].round(1)
    expected_summary_rows=700 if args.raw_metrics else 720
    assert len(seeds)==2100 and len(summary)==expected_summary_rows and len(scales)==24
    if args.preview_dir:
        args.preview_dir.mkdir(parents=True,exist_ok=True)
    audits=[]
    with PdfPages(combined_file) as combined, PdfPages(main_file) as main, PdfPages(supp_file) as supp:
        for dataset,number in DATASETS:
            fig,bbox,audit=render_page(dataset,number,seeds,summary,scales,args.raw_metrics)
            audits.append(audit)
            combined.savefig(fig,facecolor='white')
            (main if number==1 else supp).savefig(fig,facecolor='white')
            if number==1:
                # Crop the existing composite layout in page coordinates; no reflow.
                # Keep this Illustrator placement PDF narrower than 174 mm.
                pad_pixels=1.5*fig.dpi/72
                crop=Bbox.from_extents(bbox.x0-pad_pixels,bbox.y0-pad_pixels,
                    bbox.x1+pad_pixels,bbox.y1+pad_pixels).transformed(fig.dpi_scale_trans.inverted())
                assert crop.width*25.4<=174.0
                fig.savefig(crop_file,format='pdf',facecolor='white',bbox_inches=crop,pad_inches=0)
            if args.preview_dir:
                fig.savefig(args.preview_dir/f'Dataset{number}_A4_preview.png',dpi=200,facecolor='white')
            plt.close(fig)
    pd.concat(audits,ignore_index=True).to_csv(audit_file,index=False,float_format='%.17g')
    print('Saved',combined_file,main_file,supp_file,crop_file,audit_file,sep='\n')
    print('All six pages have shared within-dataset metric scales and no P-value labels.')

if __name__=='__main__':
    main()
