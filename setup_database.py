#!/usr/bin/env python3
"""
Database setup script for the Digital Exoplanet Hunting Observatory.
Creates SQLite database and sets up schema.
"""

import os
import sys
import logging
from database import ObservatoryDatabase, DatabaseConfig

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)




def main():
    """Main setup function."""
    print("=" * 60)
    print("Digital Exoplanet Hunting Observatory - Database Setup")
    print("=" * 60)
    
    # Load or create configuration
    config = DatabaseConfig()
    
    # Check if config file exists
    config_file = os.path.join(os.path.dirname(__file__), 'database_config.json')
    if os.path.exists(config_file):
        print(f"Loading configuration from {config_file}")
        config.from_file(config_file)
    else:
        print("No configuration file found. Using defaults and environment variables.")
        print("You can create database_config.json to customize the database path")
    
    print(f"Database configuration:")
    print(f"  Database Path: {config.config['database']}")
    
    # Confirm setup
    response = input("Proceed with database setup? [y/N]: ")
    if response.lower() != 'y':
        print("Setup cancelled.")
        return
    
    try:
        # Step 1: Initialize database
        print("\n1. Initializing SQLite database...")
        db = ObservatoryDatabase(config)
        if not db.initialize_connection():
            print("Failed to initialize database connection. Exiting.")
            return
        
        # Step 2: Create database schema
        print("\n2. Creating database schema...")
        if not db.create_schema():
            print("Failed to create database schema. Exiting.")
            return
        
        # Step 3: Test connection
        print("\n3. Testing database connection...")
        if not db.test_connection():
            print("Database connection test failed. Exiting.")
            return
        
        print("\n" + "=" * 60)
        print("DATABASE SETUP COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"Database Path: {config.config['database']}")
        print("Database Type: SQLite")
        print("Tables: Created")
        print("Connection: Tested")
        print("\nYour observatory database is ready to use!")
        
        # Show example usage
        print("\nExample usage:")
        print("from database import ObservatoryDatabase")
        print("db = ObservatoryDatabase()")
        print("db.initialize_connection()")
        print("stats = db.get_session_statistics()")
        
        db.close()
        
    except KeyboardInterrupt:
        print("\nSetup interrupted by user.")
    except Exception as e:
        logger.error(f"Setup failed: {e}")
        print(f"Setup failed: {e}")


if __name__ == "__main__":
    main()