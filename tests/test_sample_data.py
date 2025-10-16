#!/usr/bin/env python3
"""
Sample data generator for ProcessTracker tests
Provides realistic test data and scenarios
"""

from datetime import datetime, timedelta
import random


class TestDataGenerator:
    """Generate test data for ProcessTracker tests"""
    
    @staticmethod
    def sample_srx_accessions():
        """Generate sample SRX accession numbers"""
        return [
            "SRX21843330",
            "SRX21843332", 
            "SRX18482691",
            "SRX19980238",
            "SRX20282302"
        ]
    
    @staticmethod
    def sample_organisms():
        """Generate sample organism names"""
        return ["human", "mouse"]
    
    @staticmethod
    def sample_process_types():
        """Generate sample process types"""
        return ["scRecounter", "postprocess", "QC", "analysis"]
    
    @staticmethod
    def generate_experiment_id(prefix="TEST"):
        """Generate unique experiment ID"""
        timestamp = int(datetime.now().timestamp())
        return f"{prefix}_{timestamp}_{random.randint(1000, 9999)}"
    
    @staticmethod
    def generate_process_data(count=5):
        """Generate multiple process data records"""
        data = []
        srx_list = TestDataGenerator.sample_srx_accessions()
        organisms = TestDataGenerator.sample_organisms()
        process_types = TestDataGenerator.sample_process_types()
        
        for i in range(count):
            data.append({
                'experiment_id': TestDataGenerator.generate_experiment_id(),
                'srx_accession': random.choice(srx_list),
                'organism': random.choice(organisms),
                'process_type': random.choice(process_types),
                'status': random.choice([0, 1, None]),
                'start_datetime': datetime.now() - timedelta(hours=random.randint(1, 48)),
                'finish_datetime': datetime.now() - timedelta(hours=random.randint(0, 24)) if random.choice([True, False]) else None,
                'path': f"/path/to/results/{random.randint(1000, 9999)}",
                'error_message': "Sample error message" if random.choice([True, False, False]) else None
            })
        
        return data
    
    @staticmethod
    def generate_sample_lifecycle_data():
        """Generate data representing a complete process lifecycle"""
        base_data = {
            'experiment_id': TestDataGenerator.generate_experiment_id("LIFECYCLE"),
            'srx_accession': "SRX21843330",
            'organism': "human",
            'process_type': "scRecounter"
        }
        
        # Process start
        start_data = {
            **base_data,
            'status': None,
            'start_datetime': datetime.now() - timedelta(hours=2),
            'finish_datetime': None,
            'path': None,
            'error_message': None
        }
        
        # Process completion (successful)
        finish_data = {
            **base_data,
            'status': 0,
            'start_datetime': datetime.now() - timedelta(hours=2),
            'finish_datetime': datetime.now() - timedelta(minutes=30),
            'path': "/path/to/successful/results",
            'error_message': None
        }
        
        # Process completion (failed)
        failed_data = {
            **base_data,
            'experiment_id': TestDataGenerator.generate_experiment_id("LIFECYCLE_FAIL"),
            'status': 1,
            'start_datetime': datetime.now() - timedelta(hours=1),
            'finish_datetime': datetime.now() - timedelta(minutes=15),
            'path': None,
            'error_message': "Process failed due to insufficient resources"
        }
        
        return {
            'start': start_data,
            'success': finish_data,
            'failed': failed_data
        }


def main():
    """Generate and display sample test data"""
    print("Sample Test Data for ProcessTracker")
    print("=" * 50)
    
    generator = TestDataGenerator()
    
    print("\nSample SRX Accessions:")
    for srx in generator.sample_srx_accessions():
        print(f"  - {srx}")
    
    print("\nSample Process Data:")
    for i, data in enumerate(generator.generate_process_data(3), 1):
        print(f"  {i}. Experiment: {data['experiment_id']}")
        print(f"     SRX: {data['srx_accession']}, Organism: {data['organism']}")
        print(f"     Type: {data['process_type']}, Status: {data['status']}")
    
    print("\nLifecycle Test Data:")
    lifecycle = generator.generate_sample_lifecycle_data()
    for phase, data in lifecycle.items():
        print(f"  {phase.upper()}:")
        print(f"    Experiment: {data['experiment_id']}")
        print(f"    Status: {data['status']}, Error: {data['error_message']}")


if __name__ == "__main__":
    main()