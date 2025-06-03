
from Bio import Phylo
from Bio.Phylo.Newick import Clade
import matplotlib.pyplot as plt


def load_tree(file_path):
    trees = list(Phylo.parse(file_path, "newick"))
    return trees[0]  # Ambil pohon pertama

def save_tree(tree, file_path):
    Phylo.write(tree, file_path, "newick")

def add_new_branch(tree, parent_species_name, new_species_name, branch_length=1.0):
    """
    Sisipkan cabang baru di bawah parent_species_name.
    Jika parent adalah terminal node, buat node internal baru.
    """
    for clade in tree.find_clades():
        if clade.name == parent_species_name:
            # Jika parent adalah terminal (daun)
            if clade.is_terminal():
                # Simpan nama dan branch_length lama
                old_name = clade.name
                old_branch = clade.branch_length

                # Buat clade internal baru
                new_internal = Clade(branch_length=old_branch)
                new_internal.clades.append(Clade(name=old_name, branch_length=branch_length))
                new_internal.clades.append(Clade(name=new_species_name, branch_length=branch_length))
                clade.name = None
                clade.clades = new_internal.clades
                clade.branch_length = new_internal.branch_length
                return tree
            else:
                # Parent adalah internal node, langsung tambahkan anak baru
                clade.clades.append(Clade(name=new_species_name, branch_length=branch_length))
                return tree

    raise ValueError(f"Parent species {parent_species_name} not found in tree.")

def draw_tree_on_axes(tree, ax, highlight_species_name=None):

    
    def get_label_color(label): # This receives a Clade object
      
   
        if label == highlight_species_name:
            return "red"
        return "black"

    if not tree: # Handle case ketika tree object = None
        
        if ax:
            ax.clear()
            ax.text(0.5, 0.5, "Error: Tree data tidak tersedia", ha='center', va='center', fontsize=10, color='red')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
        return ax.figure if ax else None
    # Bersihkan axes
    ax.clear()

    # Draw the tree on the provided axes
    Phylo.draw(tree, axes=ax, label_func=lambda x: x.name if x.name else "", 
               label_colors=get_label_color, do_show=False, branch_labels=None)


    ax.set_xlabel("") # Clear default x-axis label 
    ax.set_ylabel("") # Clear default y-axis label 
    ax.set_xticks([]) # Hide x-axis ticks
    ax.set_yticks([]) # Hide y-axis ticks
    ax.spines['top'].set_visible(False) # Hide spines 
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    return ax.figure # Return figure tree
