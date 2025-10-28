import pyrosetta
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.patches as patches
from adjustText import adjust_text
from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover

def find_nearby_residues(pdb_file, target_residues, nearby_residues, distance_cutoff=5.0):
    """
    Find residues within close proximity to a list of target residues.
    """
    pyrosetta.init(silent=True)  # Initialize silently
    pose = pyrosetta.pose_from_pdb(pdb_file)
    
    nearby_residues_dict = {}
    
    for target_residue in target_residues:
        target_xyz = pose.residue(target_residue).xyz("CA")
        nearby_residues_list = []
        
        for nearby_residue in nearby_residues:
            if target_residue == nearby_residue:
                continue
            
            nearby_xyz = pose.residue(nearby_residue).xyz("CA")
            distance = target_xyz.distance(nearby_xyz)
            
            if distance < distance_cutoff:
                nearby_residues_list.append((nearby_residue, distance))
        
        nearby_residues_dict[target_residue] = nearby_residues_list
    
    return nearby_residues_dict

def plot_results(result_df, motif, glycan_positions, removed_positions, 
                percentage_cutoff, ptm_cutoff, construct):
    """
    Plot the analysis results.
    """
    font_size_ax = 12
    fig, ax = plt.subplots(figsize=[10,5])
    
    x = result_df['PTMPredictionMetric']
    y = result_df['d_total_score']
    marker_size = 33

    # Define color mapping
    color_map = {
        'Potential glycan sites': 'black',
        'Introduced glycan motifs': 'black',
        'Wild-type glycosylation': 'red',
        'Too close to wt glycan': 'grey'
    }

    # Set up plotting based on glycan positions
    if not glycan_positions:
        result_df['color'] = 'Potential glycan sites'
        scatter_plot = sns.scatterplot(x=x, y=y, hue=result_df['color'],
            palette={'Potential glycan sites': color_map['Potential glycan sites']}, 
            linewidth=0, legend=True)
    else:
        result_df['color'] = 'Introduced glycan motifs'
        result_df.loc[result_df['glycan_pos'].isin(glycan_positions), 'color'] = 'Wild-type glycosylation'
        result_df.loc[result_df['glycan_pos'].isin(removed_positions), 'color'] = 'Too close to wt glycan'
        scatter_plot = sns.scatterplot(x=x, y=y, hue=result_df['color'],
            palette=color_map, linewidth=0, legend=True)

    # Add and adjust labels
    texts = []
    for i, row in result_df.iterrows():
        texts.append(plt.text(row['PTMPredictionMetric'], row['d_total_score'], 
                            str(int(row['glycan_pos'])), 
                            fontsize=8,
                            color=color_map[row['color']]))
    
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
    
    # Set up plot formatting
    plt.title(f'Results of {construct} glycan masking\n', fontsize=16, fontweight='normal', color='black')
    plt.text(0.5, 1.025, f'{motif} motif', horizontalalignment='center', verticalalignment='center', 
        transform=plt.gca().transAxes, fontsize=14)
    plt.xlabel('PTMPredictionMetric', fontsize=font_size_ax)
    plt.ylabel('dtotal_score [REU]', fontsize=font_size_ax)
    
    # Configure axis formatting
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', length=4)
    ax.tick_params(which='major', length=7)
    
    plt.yticks(fontsize=font_size_ax)
    plt.xticks(rotation=45, fontsize=font_size_ax)
    
    # Add reference lines and shading
    total_score_cutoff = result_df['d_total_score'].quantile(percentage_cutoff/100)
    xmax, xmin = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    ax.hlines(total_score_cutoff, xmin, xmax, color='black', zorder=0, alpha=0.8, linestyle='--')
    ax.vlines(ptm_cutoff, ymin, ymax, color='black', zorder=0, alpha=0.8, linestyle='--')
    
    rectangle = patches.Rectangle((ptm_cutoff, total_score_cutoff), xmax-ptm_cutoff, ymin-total_score_cutoff, 
        linewidth=1, edgecolor='none', facecolor='lightgray', label='Rectangle', zorder=0)
    plt.gca().add_patch(rectangle)
    
    plt.show()

def identify_wild_type_glycans(pdb_file, debug=False):
    """
    Identify wild-type glycans in the PDB structure using PyRosetta
    
    Args:
        pdb_file (str): Path to PDB file
        debug (bool): Enable debug output
        
    Returns:
        list: List of residue positions containing wild-type glycans
    """
    if debug:
        print(f"[DEBUG] Identifying wild-type glycans in {pdb_file}")
    
    # Initialize PyRosetta silently
    pyrosetta.init(silent=True)
    
    # Load the structure
    try:
        pose = pyrosetta.pose_from_pdb(pdb_file)
        if debug:
            print(f"[DEBUG] Loaded structure with {pose.total_residue()} residues")
    except Exception as e:
        print(f"Error loading PDB file: {e}")
        return []

    glycan_positions = []
    
    # Common glycan residue names
    glycan_residues = ['GLC', 'MAN', 'BMA', 'FUC', 'NAG', 'GAL', 'NDG', 'SIA']
    
    # Iterate through residues to find glycans
    for i in range(1, pose.total_residue() + 1):
        residue = pose.residue(i)
        res_name = residue.name()
        
        # Check if residue is a glycan
        if any(sugar in res_name for sugar in glycan_residues):
            # Get the connected residue (usually ASN)
            try:
                # Get PDB numbering
                pdb_info = pose.pdb_info()
                res_num = pdb_info.number(i)
                chain = pdb_info.chain(i)
                
                if debug:
                    print(f"[DEBUG] Found glycan {res_name} at position {res_num} chain {chain}")
                
                # Add the position if not already present
                if res_num not in glycan_positions:
                    glycan_positions.append(res_num)
                    
            except Exception as e:
                if debug:
                    print(f"[DEBUG] Error processing glycan at position {i}: {e}")
                continue
    
    if debug:
        print(f"[DEBUG] Found {len(glycan_positions)} wild-type glycans at positions: {sorted(glycan_positions)}")
    
    return sorted(glycan_positions) 