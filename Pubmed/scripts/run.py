import subprocess
import threading
import os

def run_app():
    try:
        print("Starting app.py with Uvicorn...")
        process = subprocess.Popen(['uvicorn', 'scripts.app:app', '--reload', '--port', '8080'])
        process.wait()
    except Exception as e:
        print(f"Error running app.py: {e}")

def run_ui():
    try:
        print("Starting ui.py with Streamlit...")
        ui_path = os.path.join(os.path.dirname(__file__), 'ui.py')
        process = subprocess.Popen(['streamlit', 'run', ui_path])
        process.wait()
    except Exception as e:
        print(f"Error running ui.py: {e}")

def main():
    app_process = threading.Thread(target=run_app, daemon=True)
    ui_process = threading.Thread(target=run_ui, daemon=True)

    app_process.start()
    ui_process.start()

    app_process.join()
    ui_process.join()

if __name__ == "__main__":
    main()