from dotenv import load_dotenv
load_dotenv()

from src.storage.supabase_client import db

def test_connection():
    print("Testing Supabase connection...")

    try:
        response = db.client.table("runs").select("*").limit(1).execute()
        print("Connection OK")
        print("Response:", response.data)
    except Exception as e:
        print("Connection FAILED:", e)

if __name__ == "__main__":
    test_connection()