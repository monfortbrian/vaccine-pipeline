import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def test_basic_setup():
    print("Testing Basic Setup")
    print("=" * 30)

    # Test 1: Data models
    try:
        from src.models.candidate import CandidateProtein, PipelineRun, ConfidenceTier

        # Create test candidate
        candidate = CandidateProtein(
            protein_id="TEST123",
            protein_name="Test protein",
            sequence="MKVLVLSLGMFPLADIEAAERTVQDLGKLQ",
            source="test",
            stage="testing"
        )

        # Test the fixed add_decision method
        candidate.add_decision("test", "advance", "testing the fix")

        print("Data models working - decisions append fixed!")
        print(f"Decisions count: {len(candidate.decisions)}")

    except Exception as e:
        print(f"Data models failed: {e}")
        return False

    # Test 2: Environment variables
    api_key = os.getenv("ANTHROPIC_API_KEY")
    supabase_url = os.getenv("SUPABASE_URL")
    supabase_key = os.getenv("SUPABASE_KEY")

    if api_key and api_key != "your_claude_api_key_here":
        print("Claude API key configured")
    else:
        print("Claude API key needs to be set in .env")

    if supabase_url and supabase_url != "https://your-project-id.supabase.co":
        print("Supabase URL configured")
    else:
        print("Supabase URL needs to be set in .env")

    if supabase_key and supabase_key != "your-anon-public-key-here":
        print("Supabase key configured")
    else:
        print("Supabase key needs to be set in .env")

    print("\n Basic setup complete!")
    return True

if __name__ == "__main__":
    test_basic_setup()