import sys
import os
from dotenv import load_dotenv

load_dotenv()

project_root = os.path.dirname(os.path.abspath(__file__)) 
app_directory = os.path.join(project_root, "EcoBEE_app")
sys.path.insert(0, app_directory) 


try:
    from main_app import App 
except ImportError as e:
   
    sys.exit(1) 


def start_application():
 
    try:
        app = App()  
        app.mainloop() 
    except Exception as e:
       
       
        sys.exit(1)

if __name__ == "__main__":
 
    
    start_application()