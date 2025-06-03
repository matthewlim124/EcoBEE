
import os
import json 
import random
import google.generativeai as genai
from dotenv import load_dotenv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter

#DIR
project_root =  os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(project_root, "data")
FASTA_DIR = os.path.join(DATA_DIR, "gen_fasta")
THRESHOLD_FILE = os.path.join(DATA_DIR, "threshold.json")

GENE_TYPES = ["Hsp90", "AQP", "OR"]  # Suhu, Kelembaban, Bunga

load_dotenv()
api_key_variable = os.environ.get("GEMINI_API_KEY")

with open(THRESHOLD_FILE, "r") as f:
    thresholds = json.load(f) 



try:
    
    
    from .phylogenetic import load_tree, save_tree, add_new_branch
    
except ImportError:
   
    from phylogenetic import load_tree, save_tree, add_new_branch # save_tree might be for utility, not direct UI call
    

#Helper function
def get_average_sequence(species_genes_data, gene_type):
  
    seqs = []
    if species_genes_data:
        for species_name, genes_dict in species_genes_data.items():
            if gene_type in genes_dict and hasattr(genes_dict[gene_type], 'seq'):
                seqs.append(str(genes_dict[gene_type].seq))
            


    if not seqs:
       
        return SeqRecord(Seq(""), id=f"{gene_type}_consensus_empty", description="No sequences available")

    min_len = min(len(s) for s in seqs)
    if min_len == 0:
       
        return SeqRecord(Seq(""), id=f"{gene_type}_consensus_zerolen", description="Zero length sequences")
        
    seqs = [s[:min_len] for s in seqs]

    consensus_char_list = []
    for i in range(min_len):
        bases_at_position_i = [s[i] for s in seqs]
        if bases_at_position_i:
            most_common_base = Counter(bases_at_position_i).most_common(1)[0][0]
            consensus_char_list.append(most_common_base)
        else: 
            consensus_char_list.append("N") 

    consensus_seq_str = "".join(consensus_char_list)
    return SeqRecord(Seq(consensus_seq_str), id=f"{gene_type}_consensus", description="Consensus sequence")

#Main Simulation Function 
def run_evolution_simulation(input_env_params, filo_tree_path, all_species_genes_data):
    

    current_tree = load_tree(filo_tree_path)
    if not current_tree:
        return {"status": "error", "message": f"Failed to load base tree from {filo_tree_path}."}
    if not all_species_genes_data:
         return {"status": "error", "message": "Species gene data not provided or is empty."}

  
    matching_species = find_matching_species(input_env_params) # from rule_based.py
    if matching_species:
        
        return {
            "status": "match_found",
            "matching_species": matching_species,
            "updated_tree_object": current_tree, # No change to tree
            "parent_of_evolution": None,
            "mutated_genes_list": [],
            "new_fasta_sequences": {},
            "new_species_name": None,
        }

   
    targeted_genes = get_mutation_targets(input_env_params) # from rule_based.py
    if not targeted_genes:
       
        return {
            "status": "no_targets_for_mutation",
            "updated_tree_object": current_tree, # No change to tree
            "parent_of_evolution": None,
            "mutated_genes_list": [],
            "new_fasta_sequences": {},
            "new_species_name": None,
        }
   
  
    mutated_gene_seq_records = {} 
    new_fasta_strings = {}     
    for gene_type in targeted_genes:
        parent_seq_for_mutation = get_average_sequence(all_species_genes_data, gene_type)
        if not parent_seq_for_mutation.seq: 
           
            continue
        
        mutated_record = mutate_sequence(parent_seq_for_mutation) # from mutation.py
        mutated_gene_seq_records[gene_type] = mutated_record
        new_fasta_strings[gene_type] = f">{mutated_record.id}\n{str(mutated_record.seq)}"
       

    if not mutated_gene_seq_records: # Check if any gene was actually mutated
     
        return {
            "status": "no_genes_mutated",
            "mutated_genes_list": targeted_genes, # Still report what was targeted
            "updated_tree_object": current_tree,
            "parent_of_evolution": None,
            "new_fasta_sequences": {},
            "new_species_name": None,
        }


    fitness_scores_for_each_mutated_gene = {}
    for gene_type, m_seq_rec in mutated_gene_seq_records.items():
        scores = compute_fitness(m_seq_rec, all_species_genes_data) # from fitness.py
        fitness_scores_for_each_mutated_gene[gene_type] = scores

    species_overall_fitness = {}
    for species_name in all_species_genes_data.keys():
        current_species_total_fitness = 0
        num_genes_considered = 0
        for gene_type in mutated_gene_seq_records.keys(): # Iterate over genes that were actually mutated
            if gene_type in fitness_scores_for_each_mutated_gene and \
               species_name in fitness_scores_for_each_mutated_gene[gene_type]:
                current_species_total_fitness += fitness_scores_for_each_mutated_gene[gene_type][species_name]
                num_genes_considered += 1
        if num_genes_considered > 0:
            species_overall_fitness[species_name] = current_species_total_fitness / num_genes_considered

    parent_species_for_evolution = None
    if species_overall_fitness:
        parent_species_for_evolution = max(species_overall_fitness, key=species_overall_fitness.get)
        
    else:

        parent_species_for_evolution = list(all_species_genes_data.keys())[0] if all_species_genes_data else "UnknownParent"
       

  
    new_species_name = f"Evolved_{parent_species_for_evolution.replace('Apis_', '')}_{len(current_tree.get_terminals()) + 1}"
    
    

 
    if parent_species_for_evolution and not parent_species_for_evolution.lower().startswith("apis_"):
        
        species_part = parent_species_for_evolution.lower()
        parent_species_for_evolution = f"Apis_{species_part}"
        
    elif parent_species_for_evolution and parent_species_for_evolution.lower().startswith("apis_"):
       
        species_part = parent_species_for_evolution.split('_', 1)[1].lower()
        parent_species_for_evolution = f"Apis_{species_part}" 


    if parent_species_for_evolution:
         print(f"Final parent for new branch: {parent_species_for_evolution}")
    else: 
       
        return {"status": "error", "message": "Failed to determine a valid parent species."}


    
    
    try:
        updated_tree = add_new_branch(current_tree, parent_species_for_evolution, new_species_name, branch_length=0.1) # from tree_utils.py
    except ValueError as e: 
       
        return {
            "status": "error_updating_tree",
            "message": str(e),
            "parent_of_evolution": parent_species_for_evolution,
            "mutated_genes_list": list(mutated_gene_seq_records.keys()), 
            "new_fasta_sequences": new_fasta_strings,
            "new_species_name": new_species_name, 
            "updated_tree_object": current_tree, # Return original tree
        }
    
    
   
    return {
        "status": "evolution_simulated",
        "parent_of_evolution": parent_species_for_evolution,
        "mutated_genes_list": list(mutated_gene_seq_records.keys()), 
        "new_fasta_sequences": new_fasta_strings,
        "new_species_name": new_species_name,
        "updated_tree_object": updated_tree,
    }

#Additional Functions
def hamming_distance(seq1, seq2):
   
    if len(seq1) != len(seq2):
        raise ValueError("Sequences panjang harus sama untuk Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def compute_fitness(mutated_seq_record, species_genes):
   
    fitness_scores = {}
    for species, genes in species_genes.items():
        # genes is dict: gene_type -> SeqRecord
        try:
            dist = hamming_distance(str(mutated_seq_record.seq), str(genes[mutated_seq_record.id.split('_')[0]].seq))
            fitness = 1 - dist / len(mutated_seq_record.seq)
            fitness_scores[species] = fitness
        except Exception:
           
            continue
    return fitness_scores

def load_genes():
   
    genes = {}
    for gene_type in GENE_TYPES:
        folder = os.path.join(FASTA_DIR, gene_type)
        if not os.path.exists(folder):
            continue
        for fasta_file in os.listdir(folder):
            if not (fasta_file.endswith(".fasta") or fasta_file.endswith(".fna")):
                continue
            
            path = os.path.join(folder, fasta_file)
            
            # Misal filename: Apis_mellifera.fasta
            species_name = fasta_file.replace(".fasta", "").replace(".fna", "")
            
            record = SeqIO.read(path, "fasta")
            
            if species_name not in genes:
                genes[species_name] = {}
            
            genes[species_name][gene_type] = record
            
    return genes



def find_matching_species(input_env):
   
    for species, thr in thresholds.items():
        temp_ok = thr["temperature"][0] <= input_env["temperature"] <= thr["temperature"][1]
        hum_ok = thr["humidity"][0] <= input_env["humidity"] <= thr["humidity"][1]
        flow_ok = thr["flowers"][0] <= input_env["flowers"] <= thr["flowers"][1]
        if temp_ok and hum_ok and flow_ok:
            return species
    return None

def get_mutation_targets(input_env):
   
    # Cari species yang paling "dekat" dengan input_env, hitung berapa parameter yang cocok
    best_species = None
    best_score = -1  # jumlah parameter cocok maksimal adalah 3

    for species, thr in thresholds.items():
        score = 0
        if thr["temperature"][0] <= input_env["temperature"] <= thr["temperature"][1]:
            score += 1
        if thr["humidity"][0] <= input_env["humidity"] <= thr["humidity"][1]:
            score += 1
        if thr["flowers"][0] <= input_env["flowers"] <= thr["flowers"][1]:
            score += 1

        if score > best_score:
            best_score = score
            best_species = species

    targets = []

   
    if best_species is not None:
        thr = thresholds[best_species]
        if not (thr["temperature"][0] <= input_env["temperature"] <= thr["temperature"][1]):
            targets.append("Hsp90")
        if not (thr["humidity"][0] <= input_env["humidity"] <= thr["humidity"][1]):
            targets.append("AQP")
        if not (thr["flowers"][0] <= input_env["flowers"] <= thr["flowers"][1]):
            targets.append("OR")

    return targets

def mutate_sequence(seq_record, mutation_rate=0.05):
   
    bases = ["A", "T", "C", "G"]
    seq_str = str(seq_record.seq)
    seq_list = list(seq_str)
    n_mutations = max(1, int(len(seq_list) * mutation_rate))

    for _ in range(n_mutations):
        idx = random.randint(0, len(seq_list) - 1)
        current_base = seq_list[idx]
        possible_bases = [b for b in bases if b != current_base]
        seq_list[idx] = random.choice(possible_bases)

    mutated_seq = "".join(seq_list)
    mutated_record = SeqRecord(Seq(mutated_seq), id=seq_record.id + "_mutated", description="Mutated sequence")
    return mutated_record

def get_evolution_explanation_from_ai(prompt_data):
    if not api_key_variable:
        return "AI service not configured: API key missing."

  
    prompt = prompt_data
    
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-05-20') # Or the specific model you intend to use
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
     
        return f"Error getting explanation from AI: {e}"
    

def generate_ai_prompt(prompt_data, received_params):
    env = received_params
    mut_genes = prompt_data.get("mutated_genes", [])
    
    prompt_lines = [
        "Anda adalah seorang ahli biologi evolusioner yang menjelaskan adaptasi lebah madu.",
        "Sebuah simulasi evolusi telah dijalankan dengan kondisi lingkungan sebagai berikut:",
        f"- Suhu: {env.get('temperature', 'N/A')} °C",
        f"- Kelembapan: {env.get('humidity', 'N/A')} %",
        f"- Ketersediaan Flora: {env.get('flowers', 'N/A')} unit",
        #f"- Tingkat Predator: {env.get('predator_level', 'N/A')} indeks",
        #f"- Tingkat Polusi: {env.get('pollution_level', 'N/A')} indeks",
        "\n",
        f"Sebagai respons terhadap kondisi ini, spesies lebah '{prompt_data.get('parent_species', 'N/A')}' ",
        f"disimulasikan telah berevolusi menjadi garis keturunan baru bernama '{prompt_data.get('evolved_species_name', 'N/A')}'."
    ]

    if mut_genes:
        prompt_lines.append(f"Perubahan genetik utama yang diamati adalah mutasi pada gen berikut: {', '.join(mut_genes)}.")
        
    else:
        prompt_lines.append("Tidak ada gen spesifik yang ditargetkan untuk mutasi dalam skenario ini.")

    prompt_lines.extend([
        "\n",
        "Mohon berikan penjelasan ilmiah yang ringkas, mudah dipahami, dan langsung tanpa kata-kata yang menunjukkan Anda yang menjawab (dalam Bahasa Indonesia):",
        "1. Mengapa adaptasi evolusioner ini (melibatkan mutasi pada gen-gen tersebut, jika ada) kemungkinan besar terjadi sebagai respons terhadap kondisi lingkungan spesifik yang diberikan?",
        "2. Apa potensi keuntungan adaptif dari perubahan genetik ini bagi kelangsungan hidup spesies baru tersebut di lingkungan tersebut?",
        "Jawab kedua pertanyaan dalam 3 Kalimat !!"
    ])
    return "\n".join(prompt_lines)