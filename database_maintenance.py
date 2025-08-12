#!/usr/bin/env python3
"""
Database maintenance script for the Digital Exoplanet Hunting Observatory.
Performs routine maintenance tasks, data validation, and cleanup operations.
"""

import argparse
import logging
import sys
import schedule
import time
from datetime import datetime
from database import ObservatoryDatabase

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def perform_full_maintenance(database: ObservatoryDatabase):
    """
    Perform comprehensive database maintenance.
    
    Args:
        database: Database instance
    """
    logger.info("=" * 60)
    logger.info("STARTING FULL DATABASE MAINTENANCE")
    logger.info("=" * 60)
    
    try:
        # 1. Get current status
        logger.info("1. Checking database status...")
        status = database.get_maintenance_status()
        if not status.empty:
            logger.info(f"Database tables status:\n{status}")
        
        # 2. Validate and promote staging data
        logger.info("2. Validating and promoting staging data...")
        staging_stats = database.validate_and_promote_staging()
        if 'error' not in staging_stats:
            logger.info(f"Staging promotion: {staging_stats}")
        else:
            logger.error(f"Staging promotion failed: {staging_stats['error']}")
        
        # 3. Clean old staging data
        logger.info("3. Cleaning up old staging data...")
        cleanup_result = database.cleanup_staging_data(days_old=7)
        logger.info(f"Staging cleanup: {cleanup_result}")
        
        # 4. Perform VACUUM/ANALYZE
        logger.info("4. Performing database maintenance (VACUUM/ANALYZE)...")
        maintenance_result = database.perform_maintenance()
        logger.info(f"Database maintenance: {maintenance_result}")
        
        # 5. Final status check
        logger.info("5. Final status check...")
        final_status = database.get_maintenance_status()
        if not final_status.empty:
            logger.info(f"Final database status:\n{final_status}")
        
        logger.info("=" * 60)
        logger.info("FULL DATABASE MAINTENANCE COMPLETED SUCCESSFULLY")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"Full maintenance failed: {e}")
        return False
    
    return True


def perform_light_maintenance(database: ObservatoryDatabase):
    """
    Perform light maintenance (staging cleanup only).
    
    Args:
        database: Database instance
    """
    logger.info("Performing light database maintenance...")
    
    try:
        # Clean staging data older than 3 days
        cleanup_result = database.cleanup_staging_data(days_old=3)
        logger.info(f"Light maintenance completed: {cleanup_result}")
        return True
        
    except Exception as e:
        logger.error(f"Light maintenance failed: {e}")
        return False


def check_database_health(database: ObservatoryDatabase):
    """
    Check database health and report issues.
    
    Args:
        database: Database instance
    """
    logger.info("Checking database health...")
    
    try:
        # Test connection
        if not database.test_connection():
            logger.error("Database connection test failed")
            return False
        
        # Get maintenance status
        status = database.get_maintenance_status()
        if status.empty:
            logger.warning("Could not retrieve maintenance status")
            return False
        
        # Check for large staging tables
        for _, row in status.iterrows():
            if 'staging' in row['table_name']:
                # Parse size to check if it's getting large
                size_str = row['size']
                if 'GB' in size_str:
                    size_val = float(size_str.split()[0])
                    if size_val > 1.0:  # > 1GB
                        logger.warning(f"Large staging table detected: {row['table_name']} ({size_str})")
        
        logger.info("Database health check completed")
        return True
        
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return False


def run_maintenance_daemon():
    """Run maintenance as a daemon with scheduling."""
    logger.info("Starting database maintenance daemon...")
    
    database = ObservatoryDatabase()
    if not database.initialize_connection_pool():
        logger.error("Failed to initialize database connection")
        return 1
    
    try:
        # Schedule maintenance tasks
        schedule.every(1).hours.do(check_database_health, database)
        schedule.every(6).hours.do(perform_light_maintenance, database)
        schedule.every().day.at("02:00").do(perform_full_maintenance, database)
        
        logger.info("Scheduled maintenance tasks:")
        logger.info("  - Health check: Every hour")
        logger.info("  - Light cleanup: Every 6 hours") 
        logger.info("  - Full maintenance: Daily at 2:00 AM")
        
        # Run scheduled tasks
        while True:
            schedule.run_pending()
            time.sleep(60)  # Check every minute
            
    except KeyboardInterrupt:
        logger.info("Maintenance daemon stopped by user")
        return 0
    except Exception as e:
        logger.error(f"Maintenance daemon error: {e}")
        return 1
    finally:
        database.close()


def main():
    """Main maintenance script entry point."""
    parser = argparse.ArgumentParser(
        description="Database maintenance for Digital Exoplanet Hunting Observatory"
    )
    
    parser.add_argument(
        "--mode", 
        choices=['full', 'light', 'health', 'daemon'],
        default='full',
        help="Maintenance mode (default: full)"
    )
    
    parser.add_argument(
        "--cleanup-days", 
        type=int, 
        default=7,
        help="Days to keep staging data (default: 7)"
    )
    
    parser.add_argument(
        "--verbose", 
        action='store_true',
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize database
    database = ObservatoryDatabase()
    if not database.initialize_connection_pool():
        logger.error("Failed to initialize database connection")
        return 1
    
    try:
        success = True
        
        if args.mode == 'full':
            success = perform_full_maintenance(database)
            
        elif args.mode == 'light':
            success = perform_light_maintenance(database)
            
        elif args.mode == 'health':
            success = check_database_health(database)
            
        elif args.mode == 'daemon':
            return run_maintenance_daemon()
        
        return 0 if success else 1
        
    except KeyboardInterrupt:
        logger.info("Maintenance interrupted by user")
        return 0
    except Exception as e:
        logger.error(f"Maintenance failed: {e}")
        return 1
    finally:
        database.close()


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)