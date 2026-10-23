"""Execute every additional computational figure analysis from the two papers.

This is a coverage runner, not a declaration of historical equivalence. Exact
author arrays, plotting choices and several normalization conventions are not
available. The report records those limits instead of fitting away differences.
All inputs are bundled; no network, training or GPU is used.
"""
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.stats import skew
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from enm import (read_ca_csv, coordinates, kirchhoff, hessian, contacts,
                 modes_and_covariance, fluctuations, scale_b_factors)
from analysis_methods import (cross_correlation, chain_matrix, principal_frame,
                              directional_variances, equilibrium_matrix)

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/extended'
KBT_300=0.00198720425864083*300  # kcal/mol, illustrative assumption, not in paper


def table(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)


def load(case):
    rows=read_ca_csv(ROOT/f'data/ca/{case}.csv')
    return rows,coordinates(rows),np.array([r['b'] for r in rows])


def gnm_analyses():
    rows,xyz,b=load('5pti_A')
    _,_,cov,_=modes_and_covariance(kirchhoff(xyz,7))
    corr=cross_correlation(cov)
    np.savez_compressed(OUT/'5pti_gnm_covariance.npz',covariance=cov,correlation=corr)
    table('5pti_correlation_pairs.csv',[dict(residue_i=rows[i]['resseq'],residue_j=rows[j]['resseq'],correlation=float(corr[i,j])) for i in range(len(rows)) for j in range(i,len(rows))])
    fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
    x=np.array([r['resseq'] for r in rows]);xx,yy=np.meshgrid(x,x)
    masked_positive=np.ma.masked_where(xx>=yy,corr)
    masked_negative=np.ma.masked_where(xx<=yy,corr)
    ax.contour(x,x,masked_positive,levels=[.21,.32,.43,.54],colors='#0a7a80',linewidths=.8)
    ax.contour(x,x,masked_negative,levels=[-.26,-.22,-.18,-.14],colors='#a43159',linewidths=.8)
    ax.plot([1,58],[1,58],color='gray',lw=.5)
    ax.invert_yaxis()
    ax.set(xlabel='5PTI residue',ylabel='5PTI residue',title='GNM Figure 4 analysis: published contour levels\nPositive lower left; negative upper right',aspect='equal')
    # Approximate NMA-region centers marked by X in the original raster. These
    # are location comparisons; the original NMA covariance array is unavailable.
    centers=[('positive',288,485),('positive',803,744),('negative',402,200),
             ('negative',760,210),('negative',764,495),('negative',800,58)]
    landmarks=[]
    for sign,px,py in centers:
        rx=(px-31)*58/(895-31);ry=(py-20)*58/(825-20)
        i=int(round(rx))-1;j=int(round(ry))-1;value=float(corr[j,i])
        compatible=value>=.21 if sign=='positive' else value<=-.14
        landmarks.append(dict(source='GNM Figure 4 X marker, earlier NMA region',
            sign=sign,pixel_x=px,pixel_y=py,residue_x_approx=rx,residue_y_approx=ry,
            compared_residue_x=i+1,compared_residue_y=j+1,GNM_correlation=value,
            within_published_GNM_outer_contour=compatible,residue_uncertainty=1))
        ax.plot(rx,ry,'x',color='black',ms=7,mew=1)
    table('gnm_figure4_NMA_region_comparison.csv',landmarks)
    save(fig,'gnm_figure4_correlations')
    mask1=(x>=16)&(x<=24);mask2=(x>=29)&(x<=35);maskh=(x>=47)&(x<=56)
    upper=np.triu_indices(maskh.sum(),1)
    claims=[dict(claim='5PTI strands 16-24 / 29-35 positively correlated',measured=float(corr[np.ix_(mask1,mask2)].mean()),criterion='mean > 0',passed=bool(corr[np.ix_(mask1,mask2)].mean()>0)),
            dict(claim='5PTI helix 47-56 internally positively correlated',measured=float(corr[np.ix_(maskh,maskh)][upper].mean()),criterion='off-diagonal mean > 0',passed=bool(corr[np.ix_(maskh,maskh)][upper].mean()>0))]
    table('gnm_text_claims.csv',claims)
    chainrows=[]
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for col,case in enumerate(('1lts_C','1omf_blank')):
        rows,xyz,b=load(case);n=len(rows);full=kirchhoff(xyz,7);path=chain_matrix(n)
        matrices={'all contacts':full,'nonbonded: full - free path':full-path,
                  'free chain':path,'literal printed Eq.8':chain_matrix(n,'literal'),
                  'ring wording':chain_matrix(n,'ring')}
        _,_,c0,_=modes_and_covariance(full);_,shared=scale_b_factors(np.diag(c0),b)
        for name,k in matrices.items():
            vals,_,c,z=modes_and_covariance(k)
            pred,own=scale_b_factors(np.diag(c),b)
            for i,r in enumerate(rows):
                chainrows.append(dict(case=case,component=name,residue=r['resseq'],zero_modes=z,
                    separately_area_calibrated_B_A2=float(pred[i]),
                    full_network_scale_B_A2=float(8*np.pi**2/3*shared*c[i,i]),
                    own_multiplier=own,full_multiplier=shared))
            style='--' if name in ('literal printed Eq.8','ring wording') else '-'
            axes[0,col].plot([r['resseq'] for r in rows],pred,style,lw=1,label=name)
            if name in ('all contacts','nonbonded: full - free path','free chain'):
                axes[1,col].plot([r['resseq'] for r in rows],8*np.pi**2/3*shared*np.diag(c),lw=1,label=name)
        axes[0,col].set(title=f'{case}: each curve area calibrated',ylabel=r'B ($\AA^2$)',xlabel='PDB residue')
        axes[1,col].set(title='Same full-network scale for each component',ylabel=r'B ($\AA^2$)',xlabel='PDB residue')
    axes[0,0].legend(fontsize=7);axes[1,0].legend(fontsize=7)
    fig.suptitle('GNM Figures 5-6: connectivity and boundary-convention audit')
    save(fig,'gnm_figures5_6_connectivity');table('gnm_chain_components.csv',chainrows)
    # Count-matching sensitivity is reported separately, not substituted for input.
    rows,xyz,b=load('1tho_A');idx=[i for i,r in enumerate(rows) if r['icode']!='A']
    sensitivity=[]
    for name,indices in [('all resolved residues',np.arange(len(rows))),('omit insertion 33A (diagnostic)',np.array(idx))]:
        _,_,c,z=modes_and_covariance(kirchhoff(xyz[indices],7));pred,m=scale_b_factors(np.diag(c),b[indices])
        sensitivity.append(dict(case='1tho_A',input_policy=name,nodes=len(indices),zeros=z,pearson_r=float(np.corrcoef(pred,b[indices])[0,1]),literal_fitted_gstar=1/m))
    table('historical_input_sensitivity.csv',sensitivity)


def force_and_spectra():
    rows,xyz,b=load('1aqb_A');metrics=[];cache={}
    cutoffs=(5.,6.,7.,10.,12.,13.,14.,15.,20.,30.,70.)
    for model,builder in [('GNM',kirchhoff),('ANM',hessian)]:
        for cutoff in cutoffs:
            matrix=builder(xyz,cutoff)
            vals,vec,cov,z=modes_and_covariance(matrix,1e-8)
            # The single extra near-zero positive mode at 7 A is important.
            for tol in (1e-8,1e-6):
                if tol==1e-8: vv,cc,zz=vals,cov,z
                else: vv,_,cc,zz=modes_and_covariance(matrix,tol)
                f=fluctuations(cc,model=='ANM');pred,m=scale_b_factors(f,b)
                gamma_kbt=(3 if model=='GNM' else 1)/m
                metrics.append(dict(model=model,cutoff_A=cutoff,zero_atol=tol,zero_modes=zz,
                    gamma_over_kBT_A_minus2=gamma_kbt,gamma_at_assumed_300K_kcal_mol_A2=KBT_300*gamma_kbt,
                    anm_threefold_alternative_at_300K=(3*KBT_300/m if model=='ANM' else KBT_300*gamma_kbt),
                    pearson_r=float(np.corrcoef(pred,b)[0,1]),smallest_retained_eigenvalue=float(vv[vv>tol].min())))
            if model=='ANM':cache[cutoff]=(vals,vec,cov,z)
    table('force_constant_cutoff_sweep.csv',metrics)
    # Native embedded Figure 3 raster (846 x 620). Coordinates picked BEFORE
    # fitting; logarithmic y ticks: gamma=1 at y374, gamma=0.1 at y526.
    # A 3-pixel uncertainty gives ~4.7% gamma uncertainty. No temperature or
    # extra factor was fitted to these points.
    pixels={'GNM':[(5,150,267),(6,214,344),(7,269,403),(10,394,479),(12,457,516)],
            'ANM':[(7,269,13),(10,394,273),(12,457,339),(13,485,359),
                   (15,538,395),(20,638,447),(30,781,487)]}
    comparisons=[]
    for model,points in pixels.items():
        for cutoff,px,py in points:
            digitized=10**((374-py)/152)
            m=next(r for r in metrics if r['model']==model and r['cutoff_A']==cutoff and r['zero_atol']==1e-6)
            standard=m['gamma_at_assumed_300K_kcal_mol_A2'];alternative=m['anm_threefold_alternative_at_300K']
            comparisons.append(dict(model=model,cutoff_A=cutoff,pixel_x=px,pixel_y=py,
                figure_gamma_approx_kcal_mol_A2=digitized,pixel_uncertainty=3,
                relative_digitization_uncertainty=10**(3/152)-1,
                standard_300K_gamma=standard,threefold_ANM_alternative_300K_gamma=alternative,
                standard_relative_error=(standard-digitized)/digitized,
                alternative_relative_error=(alternative-digitized)/digitized,
                zero_atol=1e-6,zero_modes=m['zero_modes']))
    table('anm_figure3_digitized_comparison.csv',comparisons)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for model,color in [('GNM','#087a83'),('ANM','#a03a69')]:
        m=[r for r in metrics if r['model']==model and r['zero_atol']==1e-6 and r['cutoff_A']>= (5 if model=='GNM' else 7) and r['cutoff_A']<=30]
        ax.loglog([r['cutoff_A'] for r in m],[r['gamma_at_assumed_300K_kcal_mol_A2'] for r in m],'-o',color=color,label=f'{model}, mechanical convention, assumed 300 K',ms=3)
        if model=='ANM':ax.loglog([r['cutoff_A'] for r in m],[r['anm_threefold_alternative_at_300K'] for r in m],'--',color=color,label='ANM alternative: multiply scale by 3')
        p=[r for r in comparisons if r['model']==model]
        ax.scatter([r['cutoff_A'] for r in p],[r['figure_gamma_approx_kcal_mol_A2'] for r in p],facecolors='none',edgecolors=color,s=60,label=f'Published {model}, approximate digitization')
    ax.set(xlabel='Contact cutoff (A)',ylabel='Spring constant (kcal/mol/A²)',title='ANM Figure 3: explicit normalization-convention comparison')
    ax.legend(fontsize=7);save(fig,'anm_figure3_force_constants')
    fig,axes=plt.subplots(4,1,figsize=(8,11),layout='constrained');spectral=[];frequencies=[]
    for ax,cutoff in zip(axes,(10.,13.,20.,70.)):
        vals,_,cov,z=cache[cutoff];retained=vals[vals>1e-8]
        _,scale=scale_b_factors(fluctuations(cov,True),b)
        # Article omega=sqrt(gamma*lambda) convention. Time/mass conversion is
        # unspecified; do not label as Hz or experimentally measured spectra.
        freq=np.sqrt(KBT_300/scale*retained)
        counts,bins=np.histogram(freq,bins=25)
        ax.stairs(counts,bins,fill=True,color='#187e86',alpha=.75)
        ax.set(title=f'{cutoff:g} A; {len(retained)} internal modes',ylabel='Mode count',xlabel=r'Article-style $\sqrt{\gamma\lambda}$; common assumed 300 K, mass unspecified')
        spectral.append(dict(cutoff_A=cutoff,zero_modes=z,internal_modes=len(retained),frequency_skewness=float(skew(freq)),median_article_frequency=float(np.median(freq)),minimum=float(freq.min()),maximum=float(freq.max())))
        frequencies.extend(dict(cutoff_A=cutoff,nonzero_mode=i+1,eigenvalue=float(v),article_frequency=float(f)) for i,(v,f) in enumerate(zip(retained,freq)))
    fig.suptitle('ANM Figure 4: body cutoffs 10/13/20 A; caption 70 A also tested')
    save(fig,'anm_figure4_spectra');table('anm_spectra_summary.csv',spectral);table('anm_all_frequencies.csv',frequencies)
    # Quantify the article's near-invariance assertion without selecting a cutoff.
    predicted={}
    for cutoff in (12.,13.,14.,15.):
        predicted[cutoff]=scale_b_factors(fluctuations(cache[cutoff][2],True),b)[0]
    table('anm_12_to_15_cutoff_profile_comparisons.csv',[
        dict(cutoff_i=i,cutoff_j=j,profile_pearson_r=float(np.corrcoef(predicted[i],predicted[j])[0,1]),relative_L2_difference=float(np.linalg.norm(predicted[i]-predicted[j])/np.linalg.norm(predicted[i])))
        for i in predicted for j in predicted if i<j])
    return rows,xyz,b,cache[13.]


def atom_table(name):
    rows=list(csv.DictReader((ROOT/'data/supplementary'/name).open()))
    for r in rows:
        r['resseq']=int(r['resseq'])
        for key in ('x','y','z','b','occupancy'):r[key]=float(r[key])
    return rows


def draw_backbone(ax,q,colors=None,label=None):
    ax.plot(*q.T,color='#a9acb1',lw=1,label=label)
    if colors is not None:ax.scatter(*q.T,c=colors,s=6,cmap='turbo_r')
    ax.set(xlabel='X (A)',ylabel='Y (A)',zlabel='Z (A)')
    ax.set_box_aspect(np.maximum(np.ptp(q,axis=0),1))


def draw_atoms(ax,q,color,lw=1):
    """Heavy-atom sticks using explicit distance-inferred bonds (1.1-1.9 A)."""
    d=np.linalg.norm(q[:,None]-q[None,:],axis=2)
    for i,j in zip(*np.where(np.triu((d>1.1)&(d<1.9),1))):
        ax.plot(*q[[i,j]].T,color=color,lw=lw)


def directional_and_modes(rows,xyz,b,details):
    vals,vec,cov,z=details;center,axes,singular=principal_frame(xyz)
    residue=np.array([r['resseq'] for r in rows]);q=(xyz-center)@axes
    pred,multiplier=scale_b_factors(fluctuations(cov,True),b)
    physical_cov=cov*multiplier;direction=directional_variances(physical_cov,axes)
    departures=direction-direction.mean(axis=1)[:,None]
    rms=np.sqrt(direction.sum(axis=1));degree=np.diag(kirchhoff(xyz,13))
    blocks=physical_cov.reshape(len(rows),3,len(rows),3)[np.arange(len(rows)),:,np.arange(len(rows)),:]
    block_eigs=np.linalg.eigvalsh(blocks)
    table('anm_figure7_directional_variances.csv',[dict(residue=int(residue[i]),rms_A=float(rms[i]),
        X_variance_A2=float(direction[i,0]),Y_variance_A2=float(direction[i,1]),Z_variance_A2=float(direction[i,2]),
        X_departure_A2=float(departures[i,0]),Y_departure_A2=float(departures[i,1]),Z_departure_A2=float(departures[i,2]),
        smallest_to_largest_principal_variance_ratio=float(block_eigs[i,0]/block_eigs[i,2]),contacts_13A=int(degree[i])) for i in range(len(rows))])
    positive=np.flatnonzero(vals>1e-8)
    np.savez_compressed(OUT/'1aqb_anm_covariance.npz',unit_covariance=cov,
                        principal_axes=axes,center=center,
                        eigenvalues=vals,first_12_internal_modes=vec[:,positive[:12]],
                        mean_B_calibration_multiplier=multiplier)
    (OUT/'principal_frame.json').write_text(json.dumps(dict(center=center.tolist(),axes_columns=axes.tolist(),singular_values=singular.tolist(),convention='centered SVD, decreasing singular value; deterministic right handed signs'),indent=2)+'\n')
    operator,pairs=equilibrium_matrix(xyz,13)
    residual=float(np.max(np.abs(operator@operator.T-hessian(xyz,13))))
    table('anm_equilibrium_operator_check.csv',[dict(nodes=len(rows),springs=len(pairs),hessian_BBT_max_absolute_error=residual,passed=residual<1e-12)])
    fig,axesfig=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
    for j,name in enumerate('XYZ'):axesfig[0].plot(residue,departures[:,j],lw=.9,label=name)
    axesfig[0].axhline(0,color='gray',lw=.5);axesfig[0].legend(ncol=3)
    axesfig[0].set(title='ANM Figure 7: principal-axis departures from isotropy',ylabel=r'Variance departure ($\AA^2$)')
    axesfig[1].plot(residue,rms,color='#177b80');axesfig[1].set(xlabel='PDB residue',ylabel='Total RMS displacement (A)')
    save(fig,'anm_figure7_anisotropy')
    modes=[];mode_metrics=[]
    for number,j in enumerate(positive[:12],1):
        u=vec[:,j].reshape(-1,3);u=u if u[np.unravel_index(np.abs(u).argmax(),u.shape)]>=0 else -u
        projected=u@axes;thermal=u*np.sqrt(multiplier/vals[j])
        amplitude=np.linalg.norm(u,axis=1);thermal_norm=np.linalg.norm(thermal,axis=1)
        i65=int(np.flatnonzero(residue==65)[0]);i95=int(np.flatnonzero(residue==95)[0])
        mode_metrics.append(dict(nonzero_mode=number,eigenvalue=float(vals[j]),maximum_amplitude_residue=int(residue[amplitude.argmax()]),
            residue65_unit_mode_norm=float(amplitude[i65]),residue95_unit_mode_norm=float(amplitude[i95]),
            residue95_over_65_amplitude_ratio=float(amplitude[i95]/amplitude[i65]),
            residue65_thermal_RMS_A=float(thermal_norm[i65]),residue95_thermal_RMS_A=float(thermal_norm[i95]),
            residue65_two_RMS_span_A=float(2*thermal_norm[i65]),residue95_two_RMS_span_A=float(2*thermal_norm[i95])))
        for i,r in enumerate(rows):
            modes.append(dict(nonzero_mode=number,residue=r['resseq'],eigenvalue=float(vals[j]),
                unit_X=float(projected[i,0]),unit_Y=float(projected[i,1]),unit_Z=float(projected[i,2]),
                thermal_X_A=float(thermal[i]@axes[:,0]),thermal_Y_A=float(thermal[i]@axes[:,1]),thermal_Z_A=float(thermal[i]@axes[:,2]),unit_norm=float(amplitude[i])))
    table('anm_first12_mode_vectors.csv',modes);table('anm_mode_amplitude_comparison.csv',mode_metrics)
    fig,ax=plt.subplots(figsize=(11,4),layout='constrained')
    second=[r for r in modes if r['nonzero_mode']==2]
    for name in 'XYZ':ax.plot(residue,[r['unit_'+name] for r in second],lw=.9,label=name)
    ax.legend(ncol=3);ax.set(xlabel='PDB residue',ylabel='Unit-norm mode component',title='ANM Figure 8: second internal mode, principal-axis frame')
    save(fig,'anm_figure8_global_mode')
    # Structural panels show the computational CA model, not an atomistic ribbon.
    ligand=coordinates(atom_table('1aqb_retinol_atoms.csv'));ligand=(ligand-center)@axes
    fig=plt.figure(figsize=(14,10),layout='constrained')
    a=fig.add_subplot(221,projection='3d');draw_backbone(a,q,colors=pred)
    draw_atoms(a,ligand,'black',2)
    for j,name in enumerate('XYZ'):
        v=np.zeros(3);v[j]=25;a.plot([0,v[0]],[0,v[1]],[0,v[2]],lw=1);a.text(*v,name)
    a.set_title('Figure 6 analysis: computed B-factor colors / PCA axes');a.view_init(20,35)
    # Fig1 contact network at6A, explicit central-spring displacement schematic.
    a=fig.add_subplot(222,projection='3d');draw_backbone(a,q)
    _,_,_,adj=contacts(xyz,6)
    for i,j in zip(*np.where(np.triu(adj,1))):a.plot(*q[[i,j]].T,color='#168289',alpha=.3,lw=.45)
    a.set_title('Figure 1: same CA spring-network construction, 6 A');a.view_init(20,35)
    loops=[(48,52,'blue'),(63,68,'gold'),(94,99,'green'),(120,127,'cyan')]
    sidechains=atom_table('1aqb_loop_atoms.csv')
    for pane,elev,azim in [(223,20,35),(224,85,35)]:
        a=fig.add_subplot(pane,projection='3d');draw_backbone(a,q)
        for lo,hi,color in loops:
            mask=(residue>=lo)&(residue<=hi);a.plot(*q[mask].T,color=color,lw=4,label=f'{lo}-{hi}')
        for lo,hi,color in loops:
            if lo not in (63,94):continue
            heavy=coordinates([r for r in sidechains if lo<=r['resseq']<=hi])
            draw_atoms(a,(heavy-center)@axes,color,.7)
        draw_atoms(a,ligand,'red',2);a.view_init(elev,azim);a.set_title('Figure 9 analysis: reported loop regions / retinol')
        if pane==224:a.legend(fontsize=7,loc='upper left')
    save(fig,'anm_structural_panels')
    j=positive[1];u=vec[:,j].reshape(-1,3)*np.sqrt(multiplier/vals[j]);u=u@axes
    fig=plt.figure(figsize=(9,7),layout='constrained');a=fig.add_subplot(111,projection='3d')
    for sign,color in [(-1,'#a33362'),(0,'#aaaaaa'),(1,'#147a87')]:
        a.plot(*(q+sign*u).T,color=color,lw=1,label=f'{sign:+d} thermal RMS mode amplitude')
    a.legend(fontsize=8);a.set_title('Figure 10a analysis: second internal mode; calibrated thermal scale')
    a.view_init(20,35);save(fig,'anm_figure10a_conformations')
    # Export actual 3D coordinates, both equilibrium and +/- one modal RMS.
    table('anm_figure10a_conformations.csv',[dict(residue=int(residue[i]),conformation=sign,X_A=float((q+sign*u)[i,0]),Y_A=float((q+sign*u)[i,1]),Z_A=float((q+sign*u)[i,2])) for sign in (-1,0,1) for i in range(len(rows))])
    # Fig2 is a geometric definition: export an actual contact pair displacement.
    i,j=pairs[len(pairs)//2];new=xyz+vec[:,positive[1]].reshape(-1,3)*np.sqrt(multiplier/vals[positive[1]])
    table('anm_figure2_contact_displacement.csv',[dict(residue=int(residue[k]),state=state,x=float(position[k,0]),y=float(position[k,1]),z=float(position[k,2])) for k in (i,j) for state,position in [('equilibrium',xyz),('plus_mode2_RMS',new)]])
    fig=plt.figure(figsize=(7,6),layout='constrained');a=fig.add_subplot(111,projection='3d')
    for state,position,color in [('Equilibrium',xyz,'gray'),('+ mode 2 RMS',new,'#147a87')]:
        a.plot(*position[[i,j]].T,'o-',color=color,label=state)
    for k in (i,j):
        delta=new[k]-xyz[k];a.quiver(*xyz[k],*delta,color='#a33362',arrow_length_ratio=.2)
        a.text(*xyz[k],str(residue[k]))
    a.legend(fontsize=8);a.set(xlabel='x (A)',ylabel='y (A)',zlabel='z (A)',title='Figure 2: equilibrium contact and modal displacements\nAn actual 13 A spring; arrows at calibrated scale')
    save(fig,'anm_figure2_displacements')
    anisotropy=np.linalg.norm(departures,axis=1)
    summary=dict(rms_max_residue=int(residue[rms.argmax()]),anisotropy_departure_max_residue=int(residue[anisotropy.argmax()]),
        anisotropy_RMS_pearson=float(np.corrcoef(anisotropy,rms)[0,1]),
        contacts_RMS_pearson=float(np.corrcoef(degree,rms)[0,1]),
        mode2_residue95_over65=mode_metrics[1]['residue95_over_65_amplitude_ratio'],published_figure10a_95_over65=8.7/6.8,
        mode2_largest_residue=mode_metrics[1]['maximum_amplitude_residue'],
        note='Published absolute amplitudes 6.8/8.7 A have an unreported scaling convention. No fit to these values was performed.')
    (OUT/'anm_claim_comparisons.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def complex_context():
    rows=atom_table('1qab_protein_atoms.csv');ca=[r for r in rows if r['atom']=='CA']
    rbp=[r for r in rows if r['chain'] in ('E','F')];ttr=[r for r in rows if r['chain'] in 'ABCD']
    other=coordinates(ttr);interface=[]
    # This 5A heavy-atom criterion is our explicit structural summary, not a
    # threshold or numerical benchmark reported in the 2001 article.
    for chain in 'EF':
        for residue in sorted({r['resseq'] for r in rbp if r['chain']==chain}):
            atoms=[r for r in rbp if r['chain']==chain and r['resseq']==residue]
            xyz=coordinates(atoms);d=np.linalg.norm(xyz[:,None,:]-other[None,:,:],axis=2)
            i,j=np.unravel_index(d.argmin(),d.shape)
            interface.append(dict(rbp_chain=chain,rbp_residue=residue,rbp_resname=atoms[0]['resname'],
                nearest_ttr_chain=ttr[j]['chain'],nearest_ttr_residue=ttr[j]['resseq'],
                nearest_heavy_atom_distance_A=float(d[i,j]),within_5A=bool(d[i,j]<=5)))
    table('1qab_RBP_TTR_heavy_atom_interface.csv',interface)
    fig=plt.figure(figsize=(9,7),layout='constrained');ax=fig.add_subplot(111,projection='3d')
    for chain in 'ABCDEF':
        q=coordinates([r for r in ca if r['chain']==chain])
        ax.plot(*q.T,color='#4e995a' if chain in 'EF' else '#9b63a6',lw=1,label=f'{"RBP" if chain in "EF" else "TTR"} {chain}')
    ax.legend(fontsize=7);ax.set_title('Figure 10b structural context: experimental human 1QAB complex\nSeparate from pig 1AQB dynamics; no binding-affinity prediction')
    save(fig,'anm_figure10b_experimental_complex')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    gnm_analyses();print('GNM correlation, connectivity and input sensitivity complete',flush=True)
    rows,xyz,b,details=force_and_spectra();print('Cutoff, spring-scale, digitization and spectra complete',flush=True)
    summary=directional_and_modes(rows,xyz,b,details);complex_context()
    manifest=dict(run_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,platform=platform.platform(),command=' '.join(sys.argv),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'enm.py',ROOT/'analysis_methods.py',ROOT/'reproduce_figures.py']},
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'data').rglob('*.csv'))})
    (OUT/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print('All additional analysis stages executed; consult figure coverage for equivalence status.')


if __name__=='__main__':main()
