# ecobee_app/main_app.py

import customtkinter as ctk
from PIL import Image # Required for CTkImage
import os # For path joining
import matplotlib as plt
import threading
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use('TkAgg') 

try:
    import simulation  
    import phylogenetic 
except ImportError as e:
   
    

    simulation = None


# Define paths untuk assets #
project_root_for_assets = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
assets_base_path = os.path.join(project_root_for_assets, "EcoBEE_app","assets", "images")
data_base_path = os.path.join(project_root_for_assets, "EcoBEE_app","data")
filo_tree_path = os.path.join(data_base_path, "BEE_prunedtree_APIS.nwk")
APP_LOGO_PATH = os.path.join(assets_base_path, "app_logo.png")
DASHBOARD_ICON_PATH = os.path.join(assets_base_path, "dashboard_icon.png")

KELUAR_ICON_PATH = os.path.join(assets_base_path, "keluar_icon.png")
DASHBOARD_ICON_ACTIVE_PATH = os.path.join(assets_base_path, "dashboard_activated_icon.png")



class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("EcoBee Dashboard")
        self.geometry("1050x850")
        
        # Matplotlib Canvas
        self.tree_figure = None
        self.tree_canvas = None
        self.filogeni_plot_frame = None # Will hold the canvas

        # Global Theme & Font Settings 
        ctk.set_appearance_mode("Light")
        self.poppins_font = "Poppins"
        self.press2p_font = "Press Start 2P"

        self.font_regular = ctk.CTkFont(family=self.poppins_font, size=16)
        self.font_regular_smaller =  ctk.CTkFont(family=self.poppins_font, size=13)
        self.font_app_name = ctk.CTkFont(family=self.press2p_font, size=16, weight="bold") # Font for "EcoBee"
        self.font_bold = ctk.CTkFont(family=self.poppins_font, size=16, weight="bold")
        self.font_logo_text = ctk.CTkFont(family=self.poppins_font, size=20, weight="bold")
        self.font_sidebar_buttons = ctk.CTkFont(family=self.poppins_font, size=17)
     
        self.font_section_titles = ctk.CTkFont(family=self.press2p_font, size=16, weight="bold")

        self.sidebar_bg_color = "#E6C562"
        self.main_bg_color = "#FAFAFA"

        self.app_name_text_color = "#000000"
        self.active_button_fg_color = "#FFFFFF"
        self.active_button_text_color = "#000000"
        self.inactive_button_text_color = "#FFFFFF"
        self.button_color_perubahan = "#ADD8E6"
        self.text_color_main = "#FFFFFF"
        
        #Initial Data 
        self.initial_species_genes_data = None
        if simulation:
            try:
                
                self.initial_species_genes_data = simulation.load_genes()
                
            except Exception as e:
                print(f"[App] Error loading initial species gene data: {e}")
               
        else:
            print("[App] Error: gen_loader module not imported. Simulation features will be affected.")

        self.filo_tree_path = filo_tree_path

        self.configure(fg_color=self.main_bg_color)

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.sidebar_frame = ctk.CTkFrame(self, width=250, corner_radius=0, fg_color=self.sidebar_bg_color)
        
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_propagate(False) #buat frame mengikuti width
        self.sidebar_frame.grid_columnconfigure(0, weight=1)
        self.sidebar_frame.grid_rowconfigure(5, weight=1)

     


        
       

        #logo config
        self.logo_image = self.load_icon(APP_LOGO_PATH, (120, 100))
        self.logo_display_label = ctk.CTkLabel(self.sidebar_frame, text="")
        if self.logo_image:
            
            self.logo_display_label.configure(image=self.logo_image)
        else:
            self.logo_placeholder = ctk.CTkLabel(self.sidebar_frame, text="EcoBee",
                                                 font=ctk.CTkFont(family=self.poppins_font, size=30, weight="bold"),
                                                 text_color="#000000")
            self.logo_placeholder.grid(row=0, column=0, padx=20, pady=(25, 15))
        self.logo_display_label.grid(row=0, column=0, padx=20, pady=(25, 15))
        self.app_name_label = ctk.CTkLabel(self.sidebar_frame, text="EcoBee",
                                           font=self.font_app_name,
                                           text_color=self.app_name_text_color)
        self.app_name_label.grid(row=1, column=0, padx=20, pady=(0, 20)) # pady=(top_spacing_to_logo, bottom_spacing_to_buttons)
        
        nav_button_start_row = 2 # Text is now in row 1
        nav_items = [
            ("Dashboard", DASHBOARD_ICON_PATH, DASHBOARD_ICON_ACTIVE_PATH,self.show_dashboard_frame),
            
          
        ]

        self.nav_buttons = {}
        self.nav_icons_normal = {} # Store normal icons
        self.nav_icons_active = {} # Store active icons

        for i, (text, icon_path_normal, icon_path_active, command) in enumerate(nav_items):
            normal_icon_image = self.load_icon(icon_path_normal, (24, 24))
            active_icon_image = self.load_icon(icon_path_active, (24, 24))

            self.nav_icons_normal[text] = normal_icon_image
            self.nav_icons_active[text] = active_icon_image

            button = ctk.CTkButton(
                self.sidebar_frame,
                text=text,
                text_color=self.inactive_button_text_color,
                font=self.font_sidebar_buttons,
                hover=False,
                anchor="w",
                image=normal_icon_image, # Start with normal icon
                command=command,
                border_width=0 
            )
            button.grid(row=i + nav_button_start_row, column=0, padx=50, pady=7, sticky="ew")
            self.nav_buttons[text] = button

        last_nav_button_row = nav_button_start_row + len(nav_items) - 1
        self.expanding_empty_row_index = last_nav_button_row + 2
        self.exit_button_row_index = self.expanding_empty_row_index + 1

        self.sidebar_frame.grid_rowconfigure(self.expanding_empty_row_index, weight=1)

        keluar_icon_image = self.load_icon(KELUAR_ICON_PATH, (24, 24))
        self.exit_button = ctk.CTkButton(
            self.sidebar_frame,
            text="Keluar",
            font=self.font_sidebar_buttons,
            image=keluar_icon_image,
            compound="left",
            anchor="w",
            fg_color="transparent",
            text_color=self.inactive_button_text_color,
            hover=False, 
            height=50,
            corner_radius=10,
            command=self.quit,
            border_width=0 
           
        )
        self.exit_button.grid(row=self.exit_button_row_index, column=0, padx=50, pady=(10, 40), sticky="sew")

        self.main_content_outer_frame = ctk.CTkFrame(self, corner_radius=0, fg_color=self.main_bg_color)
        self.main_content_outer_frame.grid(row=0, column=1, sticky="nsew")
        self.main_content_outer_frame.grid_columnconfigure(0, weight=1)
        self.main_content_outer_frame.grid_rowconfigure(1, weight=1)

 
        
       
    
        self.content_frame_bg = ctk.CTkFrame(self.main_content_outer_frame, corner_radius=25, fg_color="transparent")
        self.content_frame_bg.grid(row=1, column=0, sticky="nsew", padx=20, pady=(0,20))
        self.content_frame_bg.grid_columnconfigure(0, weight=1)
        self.content_frame_bg.grid_rowconfigure(0, weight=1)

        self.dashboard_frame = ctk.CTkFrame(self.content_frame_bg, fg_color="transparent")
        
        

        self.current_active_button_key = "Dashboard"
        self.show_dashboard_frame()

      
   


    

    
    def load_icon(self, path, size):

        #condition untuk loading icon
        if not os.path.exists(path):
            
            return None
        try:
            return ctk.CTkImage(Image.open(path).resize(size, Image.Resampling.LANCZOS), size=size)
        except FileNotFoundError: 
           
            return None
        except Exception as e:
           
            return None

    def _update_button_styles(self, active_key):
        for key, button in self.nav_buttons.items():
            if key == active_key:
                button.configure(
                    fg_color=self.active_button_fg_color,   
                    text_color=self.active_button_text_color,
                    font=ctk.CTkFont(family=self.poppins_font, size=17, weight="bold"),
                    hover_color=self.active_button_fg_color,
                    image=self.nav_icons_active.get(key), 
                    border_width=0
                )
            else:
                button.configure(
                    fg_color='transparent',
                    text_color=self.inactive_button_text_color,
                    font=self.font_sidebar_buttons,
                    hover_color="#0000000A", 
                    image=self.nav_icons_normal.get(key) 
                )
        # Exit button
        self.exit_button.configure()


    def show_frame(self, frame_to_show, active_button_key):
        self.dashboard_frame.grid_forget()

        frame_to_show.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)
        self.current_active_button_key = active_button_key
        self._update_button_styles(active_button_key)


    def show_dashboard_frame(self):
        self.show_frame(self.dashboard_frame, "Dashboard")
        for widget in self.dashboard_frame.winfo_children(): widget.destroy()
    
        self.create_dashboard_content(self.dashboard_frame) 

    def create_dashboard_content(self, parent_frame):
      
        for widget in parent_frame.winfo_children():
            widget.destroy()

        parent_frame.grid_columnconfigure(0, weight=2)
        parent_frame.grid_columnconfigure(1, weight=1)
        parent_frame.grid_rowconfigure(0, weight=1) 
        parent_frame.grid_rowconfigure(1, weight=2) 

       

        # Filogeni Section Frame (Top-Left: Row 0, Column 0)
        filogeni_outer_frame = ctk.CTkFrame(parent_frame, corner_radius=10, fg_color="#FFFFFF")
        filogeni_outer_frame.grid(row=0, column=0, sticky="nsew", pady=(20,10))
        filogeni_outer_frame.grid_columnconfigure(0, weight=1)
        filogeni_outer_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(filogeni_outer_frame, text="Filogeni", font=self.font_section_titles).grid(row=0, column=0, pady=(10,5), padx=20, sticky="nw")
        self.filogeni_plot_frame = ctk.CTkFrame(filogeni_outer_frame, fg_color="lightgrey") 
        self.filogeni_plot_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        self.filogeni_plot_frame.grid_columnconfigure(0, weight=1) 
        self.filogeni_plot_frame.grid_rowconfigure(0, weight=1)

        # Load and display tree
        self.display_phylogenetic_tree(tree_file_path=filo_tree_path)

       

        

        # Perubahan Parameter Section (Top-Right: Row 0, Column 1) 
        parameter_frame = ctk.CTkFrame(parent_frame, corner_radius=10, fg_color="#FFFFFF")
        parameter_frame.grid(row=0, column=1, sticky="nsew", padx=(10, 0), pady=(20,10))
        parameter_frame.grid_columnconfigure(0, weight=1) #content bisa expand
        ctk.CTkLabel(parameter_frame, text="Parameter Lingkungan", font=self.font_section_titles).pack(pady=10, padx=20, anchor="nw")

        param_form_frame = ctk.CTkFrame(parameter_frame, fg_color="transparent")
        param_form_frame.pack(side="top", fill="x", expand=False, padx=15, pady=5) 
        param_form_frame.grid_columnconfigure(1, weight=1) 

        self.parameter_entries = {} 

        parameters_to_set = [
            ("Suhu (°C):", "temperature"),
            ("Kelembapan (%):", "humidity"),
            ("Ketersediaan Flora (unit):", "flowers"), 
            #("Tingkat Predator (indeks):", "predator"),
            #("Tingkat Polusi (indeks):", "polusi")
        ]

        #cuman bisa ketik angka
        validate_cmd = (self.register(self.validate_numeric_input_with_exceptions), '%P') 

        for i, (label_text, key) in enumerate(parameters_to_set):
            label = ctk.CTkLabel(param_form_frame, text=label_text, font=self.font_regular_smaller)
            label.grid(row=i, column=0, padx=(5,10), pady=7, sticky="w")
            
            entry = ctk.CTkEntry(
                param_form_frame,
                font=self.font_regular_smaller,
                width=120, 
                validate="key", 
                validatecommand=validate_cmd
            )
            entry.grid(row=i, column=1, padx=5, pady=7, sticky="ew")
            self.parameter_entries[key] = entry
        
        start_button = ctk.CTkButton(
            parameter_frame,
            text="Start Simulasi",
            text_color="#000000",
            font=self.font_regular,
            hover_color="#D8F2FF",
            fg_color="#D8F2FF",
            command=self.handle_start_simulation
        )
        start_button.pack(side="top", pady=(15,15), padx=20) 
        

        #Penjelasan Evolusi Section (Now in Row 1, Column 0, Colspan 2)
        self.evolusi_explanation_frame = ctk.CTkFrame(parent_frame, corner_radius=10, fg_color="#FFFFFF")
        self.evolusi_explanation_frame.grid(row=1, column=0, columnspan=2, sticky="new", padx=5, pady=(10,0))
        self.evolusi_explanation_frame.grid_columnconfigure(0, weight=1) 
        self.evolusi_explanation_frame.grid_rowconfigure(0, weight=0) # Title
        self.evolusi_explanation_frame.grid_rowconfigure(1, weight=0) # Details labels
        self.evolusi_explanation_frame.grid_rowconfigure(2, weight=0) # FASTA Title
        self.evolusi_explanation_frame.grid_rowconfigure(3, weight=0) # FASTA Textbox 
        self.evolusi_explanation_frame.grid_rowconfigure(4, weight=0) # AI Title
        self.evolusi_explanation_frame.grid_rowconfigure(5, weight=1) # AI Textbox 

        title_label = ctk.CTkLabel(self.evolusi_explanation_frame, text="Penjelasan Evolusi", font=self.font_section_titles)
        title_label.grid(row=0, column=0, padx=20, pady=(10,10), sticky="nw")

        #Frame detail
        details_frame = ctk.CTkFrame(self.evolusi_explanation_frame, fg_color="transparent")
        details_frame.grid(row=1, column=0, sticky="new", padx=20, pady=5)
        details_frame.grid_columnconfigure(1, weight=1) 

        self.parent_species_label = ctk.CTkLabel(details_frame, text="Induk Evolusi: -", font=self.font_regular_smaller, anchor="w")
        self.parent_species_label.grid(row=0, column=0, sticky="w", pady=2)
        self.parent_species_value = ctk.CTkLabel(details_frame, text="", font=self.font_bold, anchor="w") # Value label
        self.parent_species_value.grid(row=0, column=1, sticky="ew", pady=2, padx=(5,0))


        self.mutated_genes_label = ctk.CTkLabel(details_frame, text="Gen Termutasi: -", font=self.font_regular_smaller, anchor="w")
        self.mutated_genes_label.grid(row=1, column=0, sticky="w", pady=2)
        self.mutated_genes_value = ctk.CTkLabel(details_frame, text="", font=self.font_regular_smaller, anchor="w", wraplength=300) # Value label
        self.mutated_genes_value.grid(row=1, column=1, sticky="ew", pady=2, padx=(5,0))

        ctk.CTkLabel(self.evolusi_explanation_frame, text="Sequence FASTA Hasil Mutasi:", font=self.font_regular_smaller).grid(row=2, column=0, padx=20, pady=(10,2), sticky="nw") # font_regular_smaller if you define it
        self.fasta_sequence_textbox = ctk.CTkTextbox(self.evolusi_explanation_frame, height=60, font=self.font_regular_smaller, wrap="word", border_width=1)
        self.fasta_sequence_textbox.grid(row=3, column=0, sticky="nsew", padx=20, pady=(0,10))
        self.fasta_sequence_textbox.insert("0.0", "Sequence FASTA gen mutasi akan ditampilkan di sini...")
        self.fasta_sequence_textbox.configure(state="disabled") # Read-only

        ctk.CTkLabel(self.evolusi_explanation_frame, text="Alasan Evolusi (Jawaban dari AI bisa salah, mohon dicek ulang):", font=self.font_regular_smaller).grid(row=4, column=0, padx=20, pady=(10,2), sticky="nw")
        self.ai_explanation_textbox = ctk.CTkTextbox(self.evolusi_explanation_frame, height=100, font=self.font_regular_smaller, wrap="word", border_width=1)
        self.ai_explanation_textbox.grid(row=5, column=0, sticky="sew", padx=20, pady=(0,20))
        self.ai_explanation_textbox.insert("0.0", "Penjelasan dari AI ditampilkan di sini...")
        self.ai_explanation_textbox.configure(state="disabled") # Read-only

        #loading variable
        self._ai_loading_animation_running = False
        self._ai_loading_animation_dots = 0
        self._ai_animation_job_id = None # To store the ID from self.after()


    def display_phylogenetic_tree(self, tree_object=None, tree_file_path=None, species_to_highlight=None):
        if self.filogeni_plot_frame is None:
           
            return

        
        if self.tree_canvas:
            self.tree_canvas.get_tk_widget().destroy()
            self.tree_canvas = None
        if self.tree_figure:
            
            self.tree_figure = None

        


       
        tree_to_draw = None

        if tree_object != None:
            tree_to_draw = tree_object
        elif tree_file_path:
            tree_to_draw = phylogenetic.load_tree(tree_file_path)

        
        if tree_to_draw:
            
            self.tree_figure = Figure(figsize=(6, 8), dpi=100, facecolor='white') 
            ax = self.tree_figure.add_subplot(111)

          
            phylogenetic.draw_tree_on_axes(tree_to_draw, ax, highlight_species_name=species_to_highlight)
            
            
            self.tree_figure.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0, hspace=0)


           
            self.tree_canvas = FigureCanvasTkAgg(self.tree_figure, master=self.filogeni_plot_frame)
            self.tree_canvas.draw()
            canvas_widget = self.tree_canvas.get_tk_widget()
            canvas_widget.grid(row=0, column=0, sticky="nsew") 

        else:
          
            for widget in self.filogeni_plot_frame.winfo_children(): # Clear previous content
                widget.destroy()
            error_label = ctk.CTkLabel(self.filogeni_plot_frame, text="Failed to load tree.", font=self.font_regular)
            error_label.pack(padx=10, pady=10, expand=True)

    def validate_numeric_input_with_exceptions(self, P):
        
        if P == "":  # Allow empty string 
            return True
        if P == "-": # Allow typing a negative sign
            return True
        if P.endswith(".") and P.count('.') == 1: # Allow typing a decimal point
            
            if P[:-1] == "" or P[:-1] == "-" or P[:-1].replace('-', '', 1).isdigit():
                return True
        
       
        try:
            float(P)
            return True
        except ValueError:
           
            if P.startswith('-') and P[1:].replace('.', '', 1).isdigit():
                 return True 
            return False


    def handle_start_simulation(self):
        
        retrieved_params = {}
        all_inputs_valid = True
        error_messages = [] # To collect form validation errors

       
        for key, entry_widget in self.parameter_entries.items():
            value_str = entry_widget.get().strip()
            parameter_name = key.replace("_", " ").capitalize()

            if not value_str:
                error_messages.append(f"Input untuk '{parameter_name}' kosong.")
                all_inputs_valid = False 
                retrieved_params[key] = None 
                continue

            try:
                retrieved_params[key] = float(value_str)
               
            except ValueError:
                error_messages.append(f"Input tidak valid untuk '{parameter_name}': '{value_str}' (harus angka).")
                all_inputs_valid = False
                retrieved_params[key] = None

        if not all_inputs_valid:
            
            for msg in error_messages:
                print(f"- {msg}")
           
            self.update_evolution_explanation(error_message="Input parameter tidak valid. Mohon perbaiki.")
            self.stop_ai_loading_animation(final_text="Input parameter tidak valid.")
            return 

        # Proceed if valid inputs
       
        
        if not simulation or not hasattr(simulation, 'run_evolution_simulation'):
           
            self.update_evolution_explanation(error_message="Modul simulasi tidak termuat.")
            self.stop_ai_loading_animation(final_text="Modul simulasi tidak termuat.") 
            return
        if not self.initial_species_genes_data:
           
            self.update_evolution_explanation(error_message="Data gen spesies awal tidak termuat.")
            return
        if not os.path.exists(self.filo_tree_path):
          
            self.update_evolution_explanation(error_message=f"File pohon dasar tidak ditemukan.")
            return

        self.update_evolution_explanation(loading=True) 

        try:
            # Call func dari simulation.py
            evolution_results = simulation.run_evolution_simulation(
                retrieved_params,
                self.filo_tree_path,
                self.initial_species_genes_data
            )

      
            

            prompt_data = {
                    "environment_conditions": { 
                        "temperature": evolution_results.get("temperature"),
                        "humidity":  evolution_results.get("humidity"),
                        "flora_availability":  evolution_results.get("flowers"),
                       # "predator_level":data.get("predator"),
                        #"pollution_level": data.get("polusi")
                    },
                    "parent_species":  evolution_results.get("parent_of_evolution"),
                    "evolved_species_name":  evolution_results.get("new_species_name"),
                    "mutated_genes":  evolution_results.get("mutated_genes_list"),
                   
                }
            ai_prompt = simulation.generate_ai_prompt(prompt_data, retrieved_params)
            self.update_evolution_explanation(data=evolution_results, ai_text="Sedang meminta penjelasan AI...") # Set a placeholder for AI text
            if evolution_results and evolution_results.get("status") != "error" and evolution_results.get("updated_tree_object"):
                self.display_phylogenetic_tree(
                    tree_object=evolution_results["updated_tree_object"],
                    species_to_highlight=evolution_results.get("new_species_name") or evolution_results.get("matching_species")
                )
                self.start_ai_loading_animation()
                
                #threading agar ui tidak freeze saat await response dari api llm
                def threaded_ai_call():
                    
                    
                    import time
                    time.sleep(3) 
                    ai_explanation_text = simulation.get_evolution_explanation_from_ai(ai_prompt)
                    

                    
                    self.after(0, self.stop_ai_loading_animation, ai_explanation_text) 
                                   
                
                ai_thread = threading.Thread(target=threaded_ai_call)
                ai_thread.start()

                
            elif evolution_results and evolution_results.get("status") == "error":
                 
                 ai_explanation_text = evolution_results.get("message", "Penjelasan AI tidak diperlukan untuk hasil ini.") 
                 self.stop_ai_loading_animation(final_text=ai_explanation_text)
                 
           
            elif evolution_results.get("status") == "match_found" and evolution_results.get("updated_tree_object"):
                self.display_phylogenetic_tree(tree_object=evolution_results["updated_tree_object"], species_to_highlight=evolution_results.get("matching_species"))
                ai_explanation_text = f"Spesies {evolution_results.get('matching_species')} ditemukan cocok dengan kondisi lingkungan. Tidak ada evolusi baru yang disimulasikan."
                self.stop_ai_loading_animation(final_text=ai_explanation_text)
           


        except Exception as e:
          
            import traceback
            traceback.print_exc()
            self.update_evolution_explanation(error_message=f"Terjadi kesalahan internal saat simulasi: {e}")
            self.stop_ai_loading_animation(final_text=f"Terjadi kesalahan: {e}")

    def update_evolution_explanation(self, data=None, loading=False, error_message=None, ai_text=None):
        
        self.fasta_sequence_textbox.configure(state="normal")
        self.ai_explanation_textbox.configure(state="normal")
        self.fasta_sequence_textbox.delete("0.0", "end")
        self.ai_explanation_textbox.delete("0.0", "end")


        if not loading:
            self.ai_explanation_textbox.configure(state="normal")
            self.ai_explanation_textbox.delete("0.0", "end")

        if loading:
            self.parent_species_value.configure(text="Memproses...")
            self.mutated_genes_value.configure(text="Memproses...")
            self.fasta_sequence_textbox.insert("0.0", "Memproses...")
            if not self._ai_loading_animation_running: 
                self.ai_explanation_textbox.configure(state="normal")
                self.ai_explanation_textbox.insert("0.0", "Menyiapkan simulasi...")
                self.ai_explanation_textbox.configure(state="disabled")
        elif error_message:
            self.parent_species_value.configure(text="-")
            self.mutated_genes_value.configure(text="-")
            self.fasta_sequence_textbox.insert("0.0", f"Error: {error_message}")
            if not self._ai_loading_animation_running:
             self.ai_explanation_textbox.configure(state="normal")
             self.ai_explanation_textbox.insert("0.0", "Simulasi gagal.") 
             self.ai_explanation_textbox.configure(state="disabled")
        elif data:
            if data.get("status") == "match_found":
                self.parent_species_value.configure(text=f"{data.get('matching_species', 'N/A')} (Cocok, tidak ada evolusi)")
                self.mutated_genes_value.configure(text="-")
                self.fasta_sequence_textbox.insert("0.0", "Tidak ada sekuen baru.")
                if not self._ai_loading_animation_running:
                        self.ai_explanation_textbox.configure(state="normal")
                        self.ai_explanation_textbox.insert("0.0", f"Spesies {data.get('matching_species', '')} cocok. Tidak perlu evolusi")
                        self.ai_explanation_textbox.configure(state="disabled")
            elif data.get("status") == "evolution_simulated":
                self.parent_species_value.configure(text=data.get("parent_of_evolution", "N/A"))

                mut_genes = data.get("mutated_genes_list", [])
                self.mutated_genes_value.configure(text=", ".join(mut_genes) if mut_genes else "Tidak ada")

                fasta_text = ""
                for gene, seq_fasta in data.get("new_fasta_sequences", {}).items():
                    fasta_text += f"{seq_fasta}\n\n"
                self.fasta_sequence_textbox.insert("0.0", fasta_text.strip() if fasta_text else "Tidak ada sekuen termutasi.")


                
              
                self.ai_explanation_textbox.insert("0.0", ai_text)
            else: 
                self.parent_species_value.configure(text="-")
                self.mutated_genes_value.configure(text="-")
                self.fasta_sequence_textbox.insert("0.0", "Data tidak tersedia.")
                if not self._ai_loading_animation_running:
                    self.ai_explanation_textbox.configure(state="normal")
                    self.ai_explanation_textbox.insert("0.0", "Tidak ada hasil.")
                    self.ai_explanation_textbox.configure(state="disabled")
    

            # Re-lock textboxes
            self.fasta_sequence_textbox.configure(state="disabled")
            self.ai_explanation_textbox.configure(state="disabled")

    def _animate_ai_loading(self):
        if not self._ai_loading_animation_running:
            return 

        self._ai_loading_animation_dots = (self._ai_loading_animation_dots + 1) % 5 # Cycles 0, 1, 2, 3 dots
        loading_text = "Menunggu penjelasan dari AI" + "." * self._ai_loading_animation_dots
        
        self.ai_explanation_textbox.configure(state="normal")
        self.ai_explanation_textbox.delete("0.0", "end")
        self.ai_explanation_textbox.insert("0.0", loading_text)
        self.ai_explanation_textbox.configure(state="disabled")
        
       
        self._ai_animation_job_id = self.after(500, self._animate_ai_loading) # Update every 500ms

    def start_ai_loading_animation(self):
        if self._ai_loading_animation_running:
            return
        
      
        if self._ai_animation_job_id is not None:
            self.after_cancel(self._ai_animation_job_id)
            self._ai_animation_job_id = None
            
        self._ai_loading_animation_running = True
        self._ai_loading_animation_dots = 0 
        
        self.ai_explanation_textbox.configure(state="normal")
        self.ai_explanation_textbox.delete("0.0", "end")
        
        self.ai_explanation_textbox.insert("0.0", "Memproses penjelasan AI") 
        self.ai_explanation_textbox.configure(state="disabled")
        
        self._animate_ai_loading() 

    def stop_ai_loading_animation(self, final_text=""):
        self._ai_loading_animation_running = False
        if self._ai_animation_job_id is not None:
            self.after_cancel(self._ai_animation_job_id)
            self._ai_animation_job_id = None
        
        self.ai_explanation_textbox.configure(state="normal")
        self.ai_explanation_textbox.delete("0.0", "end")
        if final_text:
            self.ai_explanation_textbox.insert("0.0", final_text)
        else:
            self.ai_explanation_textbox.insert("0.0", "Tidak ada penjelasan AI yang diterima.") # Default if no text
        self.ai_explanation_textbox.configure(state="disabled")

   
    
    

if __name__ == "__main__":
    app = App()
    app.mainloop()