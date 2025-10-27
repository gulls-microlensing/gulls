#!/usr/bin/env python3
"""
Version bumping script for Gulls.

This script updates version numbers across the codebase and creates
release commits with proper tagging. It also auto-generates RELEASE_NOTES.md
from CHANGELOG.md entries.

Usage:
    python3 scripts/bump_version.py patch    # 2.0.0 -> 2.0.1
    python3 scripts/bump_version.py minor    # 2.0.1 -> 2.1.0  
    python3 scripts/bump_version.py major    # 2.1.0 -> 3.0.0
    python3 scripts/bump_version.py release  # Create release commit and tag
    python3 scripts/bump_version.py patch --revert  # 2.0.1 -> 2.0.0

Features:
- Updates version in src/gulls.cpp, CHANGELOG.md, and documentation/conf.py
- Auto-generates RELEASE_NOTES.md from CHANGELOG.md entries
- Prompts before replacing existing RELEASE_NOTES.md with different version
- Creates release commits and tags automatically
- Handles unstaged changes intelligently
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path
from datetime import datetime

# Files that contain version information
VERSION_FILES = [
    "src/gulls.cpp",
    "CHANGELOG.md",
    "documentation/source/conf.py",
]

def get_current_version():
    """Extract current version from gulls.cpp."""
    gulls_cpp = Path("src/gulls.cpp")
    if not gulls_cpp.exists():
        raise FileNotFoundError("src/gulls.cpp not found")
    
    content = gulls_cpp.read_text()
    match = re.search(r'gulls v(\d+\.\d+\.\d+)', content)
    if not match:
        raise ValueError("Could not find version in src/gulls.cpp")
    
    return match.group(1)

def bump_version(version, bump_type, revert=False):
    """Bump version number according to semantic versioning."""
    major, minor, patch = map(int, version.split('.'))
    
    if revert:
        if bump_type == "major":
            major = max(0, major - 1)
            minor = 0
            patch = 0
        elif bump_type == "minor":
            minor = max(0, minor - 1)
            patch = 0
        elif bump_type == "patch":
            patch = max(0, patch - 1)
        else:
            raise ValueError(f"Invalid bump type: {bump_type}")
    else:
        if bump_type == "major":
            major += 1
            minor = 0
            patch = 0
        elif bump_type == "minor":
            minor += 1
            patch = 0
        elif bump_type == "patch":
            patch += 1
        else:
            raise ValueError(f"Invalid bump type: {bump_type}")
    
    return f"{major}.{minor}.{patch}"

def update_gulls_cpp(new_version):
    """Update version in src/gulls.cpp."""
    gulls_cpp = Path("src/gulls.cpp")
    content = gulls_cpp.read_text()
    
    # Update version
    content = re.sub(
        r'printf\("\\n\\n\\n\\ngulls v\d+\.\d+\.\d+\\n"\);',
        rf'printf("\\n\\n\\n\\ngulls v{new_version}\\n");',
        content
    )
    
    gulls_cpp.write_text(content)
    print(f"Updated src/gulls.cpp to version {new_version}")

def update_changelog(new_version, bump_type):
    """Add new version entry to CHANGELOG.md."""
    changelog = Path("CHANGELOG.md")
    if not changelog.exists():
        print("Warning: CHANGELOG.md not found, skipping changelog update")
        return
    
    content = changelog.read_text()
    
    # Add new version entry after the first ## [version] line
    today = datetime.now().strftime("%Y-%m-%d")
    
    new_entry = f"""## [{new_version}] - {today}

### Added
- [Add new features here]

### Changed  
- [Add changes here]

### Fixed
- [Add bug fixes here]

### Security
- [Add security fixes here]

"""
    
    # Insert after the first ## [version] line
    lines = content.split('\n')
    insert_index = 0
    for i, line in enumerate(lines):
        if line.startswith('## [') and '] -' in line:
            insert_index = i
            break
    
    lines.insert(insert_index, new_entry)
    changelog.write_text('\n'.join(lines))
    print(f"Added {new_version} entry to CHANGELOG.md")

def update_conf_py(new_version):
    """Update version in documentation/conf.py."""
    conf_py = Path("documentation/source/conf.py")
    if not conf_py.exists():
        print("Warning: documentation/source/conf.py not found, skipping")
        return
    
    content = conf_py.read_text()
    
    # Update version and release (more specific patterns)
    content = re.sub(r'version = [\'"][0-9]+\.[0-9]+\.[0-9]+[\'"]', f'version = "{new_version}"', content)
    content = re.sub(r'release = [\'"][0-9]+\.[0-9]+\.[0-9]+[\'"]', f'release = "{new_version}"', content)
    
    conf_py.write_text(content)
    print(f"Updated documentation/conf.py to version {new_version}")

def extract_changelog_entry(version):
    """Extract the changelog entry for a specific version."""
    changelog = Path("CHANGELOG.md")
    if not changelog.exists():
        return None
    
    content = changelog.read_text()
    lines = content.split('\n')
    
    # Find the version entry
    start_index = None
    end_index = None
    
    for i, line in enumerate(lines):
        if line.startswith(f"## [{version}]"):
            start_index = i
            break
    
    if start_index is None:
        return None
    
    # Find the end of this entry (next ## [version] or end of file)
    for i in range(start_index + 1, len(lines)):
        if lines[i].startswith("## [") and "] -" in lines[i]:
            end_index = i
            break
    
    if end_index is None:
        end_index = len(lines)
    
    # Extract the entry
    entry_lines = lines[start_index:end_index]
    return '\n'.join(entry_lines).strip()

def generate_release_notes_from_changelog(version):
    """Generate RELEASE_NOTES.md from CHANGELOG.md entry."""
    changelog_entry = extract_changelog_entry(version)
    if not changelog_entry:
        print(f"Warning: No changelog entry found for version {version}")
        return False
    
    # Convert changelog format to release notes format
    lines = changelog_entry.split('\n')
    release_lines = []
    
    # Replace the header
    for line in lines:
        if line.startswith(f"## [{version}]"):
            date_match = re.search(r'\] - (\d{4}-\d{2}-\d{2})', line)
            date = date_match.group(1) if date_match else "TBD"
            release_lines.append(f"# Gulls v{version} Release Notes")
            release_lines.append("")
            release_lines.append(f"**Release Date:** {date}")
            release_lines.append("")
            # Determine release type based on semantic version components
            _, minor, patch = version.split(".")
            if patch != "0":
                release_lines.append("## Patch Release")
            elif minor != "0":
                release_lines.append("## Minor Release")
            else:
                release_lines.append("## Major Release")
            release_lines.append("")
        else:
            release_lines.append(line)
    
    # Add some standard sections if they don't exist
    content = '\n'.join(release_lines)
    if "## What's New" not in content and "## What's Included" not in content:
        release_lines.append("## What's New")
        release_lines.append("")
        release_lines.append("This release includes the following changes:")
        release_lines.append("")
        release_lines.append("## What's Included")
        release_lines.append("")
        release_lines.append("- **Source code**: Complete Gulls source with CMake build system")
        release_lines.append("- **Binaries**: Linux executables (GSL fallbacks - testing only)")
        release_lines.append("- **Documentation**: Built HTML documentation")
        release_lines.append("- **Smoke test plots**: Visual proof that the release works")
        release_lines.append("")
        release_lines.append("## Getting Started")
        release_lines.append("")
        release_lines.append("1. **Install Gulls** - See the [Installation Guide](https://gulls.readthedocs.io/en/latest/install_gulls.html)")
        release_lines.append("2. **Validate your inputs** - Use `python scripts/validate_inputs.py your_file.prm`")
        release_lines.append("3. **Run simulations** - See the [Running Guide](https://gulls.readthedocs.io/en/latest/run_simulations.html)")
        release_lines.append("")
        release_lines.append("## Full Changelog")
        release_lines.append("")
        release_lines.append(f"See [CHANGELOG.md](CHANGELOG.md) for the complete list of changes.")
        release_lines.append("")
        release_lines.append("---")
        release_lines.append("")
        release_lines.append(f"**Previous Release:** v{'.'.join(version.split('.')[:-1])}.{int(version.split('.')[-1])-1 if int(version.split('.')[-1]) > 0 else '0'}")
    
    return '\n'.join(release_lines)

def update_release_notes(new_version):
    """Update or create RELEASE_NOTES.md from CHANGELOG.md."""
    release_notes_path = Path("RELEASE_NOTES.md")
    
    # Check if RELEASE_NOTES.md exists and has the right version
    if release_notes_path.exists():
        content = release_notes_path.read_text()
        if f"v{new_version}" in content:
            print(f"RELEASE_NOTES.md already exists for version {new_version}")
            return
        
        # Check if it has a different version
        version_match = re.search(r'# Gulls v(\d+\.\d+\.\d+)', content)
        if version_match:
            existing_version = version_match.group(1)
            print(f"RELEASE_NOTES.md exists for version {existing_version}, but we're releasing {new_version}")
            response = input(f"Replace RELEASE_NOTES.md with new version {new_version}? (y/N): ")
            if response.lower() not in ['y', 'yes']:
                print("Keeping existing RELEASE_NOTES.md")
                return
    
    # Generate new release notes from changelog
    print(f"Generating RELEASE_NOTES.md for version {new_version} from CHANGELOG.md...")
    release_content = generate_release_notes_from_changelog(new_version)
    
    if release_content:
        release_notes_path.write_text(release_content)
        print(f"Created/updated RELEASE_NOTES.md for version {new_version}")
        print("You can edit RELEASE_NOTES.md to customize the release notes before creating the release")
    else:
        print("Could not generate release notes from changelog")

def create_release_commit(new_version):
    """Create a release commit and tag."""
    try:
        # Check if we're in a git repository
        subprocess.run(["git", "status"], check=True, capture_output=True)
        
        # Check for unstaged changes
        try:
            result = subprocess.run(["git", "diff", "--name-only"], capture_output=True, text=True, check=True)
            unstaged_files = result.stdout.strip().split('\n') if result.stdout.strip() else []
            
            if unstaged_files:
                print(f"Found unstaged changes in {len(unstaged_files)} files:")
                for file in unstaged_files[:5]:  # Show first 5 files
                    print(f"  - {file}")
                if len(unstaged_files) > 5:
                    print(f"  ... and {len(unstaged_files) - 5} more files")
                
                response = input("Include all unstaged changes in release commit? (y/N): ")
                if response.lower() in ['y', 'yes']:
                    subprocess.run(["git", "add", "."], check=True)
                    print("Added all changes to staging area")
                else:
                    print("Only committing staged changes")
            else:
                print("No unstaged changes found")
                subprocess.run(["git", "add", "."], check=True)
        except subprocess.CalledProcessError:
            # Fallback if git diff fails
            subprocess.run(["git", "add", "."], check=True)
        
        # Create commit (only if there are changes)
        commit_msg = f"Release version {new_version}"
        try:
            subprocess.run(["git", "commit", "-m", commit_msg], check=True)
            print(f"Created release commit for version {new_version}")
        except subprocess.CalledProcessError:
            print("No changes to commit - working tree is clean")
        
        # Push any unpushed commits first
        try:
            subprocess.run(["git", "push", "origin", "HEAD"], check=True)
            print("Pushed local commits to remote")
        except subprocess.CalledProcessError:
            print("No commits to push or push failed")
        
        # Create tag (handle existing tags)
        tag_name = f"v{new_version}"
        try:
            subprocess.run(["git", "tag", "-a", tag_name, "-m", f"Release {new_version}"], check=True)
            print(f"Created tag {tag_name}")
        except subprocess.CalledProcessError:
            print(f"Tag {tag_name} already exists!")
            response = input(f"Delete existing tag {tag_name} and create new one? (y/N): ")
            if response.lower() in ['y', 'yes']:
                # Delete local tag
                subprocess.run(["git", "tag", "-d", tag_name], check=True)
                # Delete remote tag
                try:
                    subprocess.run(["git", "push", "origin", "--delete", tag_name], check=True)
                    print(f"Deleted remote tag {tag_name}")
                except subprocess.CalledProcessError:
                    print(f"Remote tag {tag_name} doesn't exist or couldn't be deleted")
                # Create new tag
                subprocess.run(["git", "tag", "-a", tag_name, "-m", f"Release {new_version}"], check=True)
                print(f"Created new tag {tag_name}")
            else:
                print("Aborting release - tag already exists")
                return
        
        # Push the tag (this triggers the release workflow)
        subprocess.run(["git", "push", "origin", tag_name], check=True)
        
        print(f"Created release commit and tag {tag_name}")
        print(f"Pushed tag {tag_name} - release workflow should trigger automatically")
        
    except subprocess.CalledProcessError as e:
        print(f"Git operations failed: {e}")
        print("Please commit changes manually")

def main():
    parser = argparse.ArgumentParser(description="Bump Gulls version number")
    parser.add_argument("bump_type", choices=["patch", "minor", "major", "release"],
                       help="Type of version bump")
    parser.add_argument("--dry-run", action="store_true",
                       help="Show what would be changed without making changes")
    parser.add_argument("--revert", action="store_true",
                       help="Revert version bump (e.g., patch --revert: 2.0.1 -> 2.0.0)")
    
    args = parser.parse_args()
    
    try:
        current_version = get_current_version()
        print(f"Current version: {current_version}")
        
        if args.bump_type == "release":
            new_version = current_version
            print(f"Creating release for version {new_version}")
        else:
            new_version = bump_version(current_version, args.bump_type, args.revert)
            if args.revert:
                print(f"Reverted version: {new_version}")
            else:
                print(f"New version: {new_version}")
        
        if args.dry_run and not args.bump_type == "release":
            print("Dry run - no changes made")
            return
        
        # Update files
        update_gulls_cpp(new_version)
        update_conf_py(new_version)
        
        # Only generate release notes when creating a release
        if args.bump_type == "release":
            update_release_notes(new_version)
            if args.dry_run:
                print("Dry run - release notes were created without git operations")
                print("You can now edit RELEASE_NOTES.md to customize the release notes")
                return
            create_release_commit(new_version)
        
        print(f"\nVersion bump complete: {current_version} -> {new_version}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
