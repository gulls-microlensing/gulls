#!/usr/bin/env python3
"""
Version bumping script for Gulls.

This script updates version numbers across the codebase and creates
release commits with proper tagging.

Usage:
    python3 scripts/bump_version.py patch    # 2.0.0 -> 2.0.1
    python3 scripts/bump_version.py minor    # 2.0.1 -> 2.1.0  
    python3 scripts/bump_version.py major    # 2.1.0 -> 3.0.0
    python3 scripts/bump_version.py release  # Create release commit and tag
    python3 scripts/bump_version.py patch --revert  # 2.0.1 -> 2.0.0
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
        
        if args.dry_run:
            print("Dry run - no changes made")
            return
        
        # Update files
        update_gulls_cpp(new_version)
        update_conf_py(new_version)
        
        if args.bump_type == "release":
            create_release_commit(new_version)
        
        print(f"\nVersion bump complete: {current_version} -> {new_version}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
