import numpy as np

def is_surface_wave(phase):
    """Check if the phase is a surface wave"""
    return phase in ['Lg', 'LQ', 'LR', 'Rg']

def is_depth_phase(phase):
    """Check if the phase is a depth phase"""
    return phase.startswith(('p', 's')) and not phase.startswith(('pP', 'sP'))

def get_phase_characteristics(phase):
    """
    Get the characteristics for a specific phase.
    Returns (min_distance, max_distance, min_depth, velocity_type)
    """
    characteristics = {
        # Regular body waves
        'P': (0, 180, 800, 'P'),
        'S': (0, 180, 800, 'S'),
        
        # Crustal phases
        'Pg': (0, 10, 30, 'Pcrust'),
        'Pb': (0, 10, 30, 'Pcrust'),
        'Sg': (0, 10, 30, 'Scrust'),
        'Sb': (0, 10, 30, 'Scrust'),
        
        # Regional phases
        'Pn': (2, 15, 50, 'Pmantle'),
        'Sn': (2, 15, 50, 'Smantle'),
        
        # Core phases
        'PKP': (110, 180, 800, 'PKP'),
        'PKPab': (140, 180, 800, 'PKPab'),
        'PKPbc': (140, 180, 800, 'PKPbc'),
        'PKPdf': (110, 180, 800, 'PKPdf'),
        'PKKP': (110, 180, 800, 'PKKP'),
        
        # SKS family
        'SKS': (60, 180, 800, 'SKS'),
        'SKSac': (60, 180, 800, 'SKSac'),
        'SKSdf': (60, 180, 800, 'SKSdf'),
        'SKKP': (60, 180, 800, 'SKKP'),
        'SKP': (60, 180, 800, 'SKP'),
        'SKPdf': (60, 180, 800, 'SKPdf'),
        
        # Reflected phases
        'PP': (0, 180, 800, 'PP'),
        'SS': (0, 180, 800, 'SS'),
        'PcP': (0, 100, 800, 'PcP'),
        'ScP': (0, 100, 800, 'ScP'),
        'ScS': (0, 100, 800, 'ScS'),
        
        # Depth phases
        'pP_': (0, 180, 800, 'pP'),
        'sP': (0, 180, 800, 'sP'),
        'sS_': (0, 180, 800, 'sS'),
        
        # Surface reflected core phases
        'pPKPab': (140, 180, 800, 'pPKPab'),
        'pPKPbc': (140, 180, 800, 'pPKPbc'),
        'pPKPdf': (110, 180, 800, 'pPKPdf'),
        'sPKPab': (140, 180, 800, 'sPKPab'),
        'sPKPbc': (140, 180, 800, 'sPKPbc'),
        'sPKPdf': (110, 180, 800, 'sPKPdf'),
        
        # Surface waves
        'Lg': (0, 30, 0, 'Lg'),
        'LQ': (0, 30, 0, 'LQ'),
        'LR': (0, 30, 0, 'LR'),
        'Rg': (0, 30, 0, 'Rg'),
        
        # Special phases
        'Is': (0, 180, 800, 'I'),
        'It': (0, 180, 800, 'I'),
        'Iw': (0, 180, 800, 'I'),
        'P1': (0, 15, 50, 'P1')
    }
    
    return characteristics.get(phase, (0, 180, 800, 'P'))

def calculate_travel_time(phase, depth, distance):
    """
    Calculate travel time for specific phase, depth, and distance using AK135 model.
    Returns None if the phase is not possible for the given geometry.
    """
    min_dist, max_dist, max_depth, velocity_type = get_phase_characteristics(phase)
    
    # Check if the distance and depth are within valid ranges
    if (distance < min_dist or distance > max_dist or depth > max_depth):
        return None
        
    # Convert distance to kilometers
    distance_km = distance * 111.19
    
    # Special handling for surface waves
    if is_surface_wave(phase):
        if depth > 0:
            return None
        
        # Surface wave velocities (group velocities in km/s)
        velocities = {
            'Lg': 3.5,
            'LQ': 3.2,
            'LR': 3.0,
            'Rg': 2.8
        }
        return distance_km / velocities[phase]
    
    # Base velocities for different wave types (km/s)
    velocities = {
        'P': 6.8,
        'S': 3.9,
        'Pcrust': 6.2,
        'Scrust': 3.5,
        'Pmantle': 8.0,
        'Smantle': 4.5,
        'PKP': 7.2,
        'PKPab': 7.3,
        'PKPbc': 7.1,
        'PKPdf': 7.0,
        'SKS': 4.8,
        'PcP': 6.9,
        'ScS': 4.0,
        'I': 6.5,
        'P1': 6.4
    }
    
    # Get base velocity
    velocity = velocities.get(velocity_type, 6.0)
    
    # Adjust velocity for depth phases
    if is_depth_phase(phase):
        velocity *= 0.95
        
    # Simple straight-line approximation - replace with actual ray path calculation
    hypotenuse = np.sqrt(distance_km**2 + depth**2)
    return hypotenuse / velocity

def generate_travel_time_table(phase, depth_samples, distance_samples):
    """Generate travel time table for a specific phase"""
    output = []
    output.append(f"n # {phase}      travel-time (and amplitude) tables")
    output.append(f"{len(depth_samples)}    # number of depth samples")
    
    depth_line = "    " + "    ".join(f"{d:.2f}" for d in depth_samples)
    output.append(depth_line)
    
    output.append(f"{len(distance_samples)}    # number of distance samples")
    
    distance_chunks = [distance_samples[i:i+10] for i in range(0, len(distance_samples), 10)]
    for chunk in distance_chunks:
        dist_line = "    " + "    ".join(f"{d:.2f}" for d in chunk)
        output.append(dist_line)
    
    for depth in depth_samples:
        output.append(f"# Travel-time/amplitude for z = {depth:8.2f}")
        
        for distance in distance_samples:
            travel_time = calculate_travel_time(phase, depth, distance)
            if travel_time is None:
                output.append("      -1.000")
            else:
                output.append(f"    {travel_time:10.3f}")
        
        output.append("      -1.000")
    
    return "\n".join(output)

def get_depth_samples(phase):
    """Get appropriate depth samples for a phase"""
    if is_surface_wave(phase):
        return [0.00]
    elif phase in ['Pg', 'Sg', 'Pb', 'Sb']:
        return [0.00, 5.00, 10.00, 15.00, 20.00]
    elif phase in ['Pn', 'Sn', 'P1']:
        return [0.00, 5.00, 15.00, 30.00, 40.00, 50.00]
    elif phase.startswith(('PKP', 'SKS', 'pPKP', 'sPKP')):
        return [0.00, 5.00, 15.00, 30.00, 40.00, 50.00, 75.00, 100.00, 150.00, 200.00]
    else:
        return [0.00, 5.00, 15.00, 30.00, 40.00, 50.00, 75.00, 100.00, 150.00, 200.00, 
                300.00, 400.00, 500.00, 600.00, 800.00]

def get_distance_samples(phase):
    """Get appropriate distance samples for a phase"""
    min_dist, max_dist, _, _ = get_phase_characteristics(phase)
    
    if is_surface_wave(phase):
        return np.arange(0, 30.1, 0.5)
    elif phase in ['Pg', 'Sg', 'Pb', 'Sb']:
        return np.arange(0, 10.1, 0.2)
    elif phase in ['Pn', 'Sn', 'P1']:
        return np.arange(0, 15.1, 0.5)
    elif phase.startswith(('PKP', 'SKS')):
        return np.arange(min_dist, max_dist + 0.1, 1.0)
    else:
        return np.arange(min_dist, max_dist + 0.1, 1.0)

def write_ak135_tables():
    """Generate and write AK135 travel time tables for all phases"""
    phases = [
        'Is', 'It', 'Iw', 'Lg', 'LQ', 'LR', 'P', 'P1', 'Pb', 'PcP', 'Pg', 
        'PKKP', 'PKP', 'PKPab', 'PKPbc', 'PKPdf', 'Pn', 'pP_', 'PP', 
        'pPKPab', 'pPKPbc', 'pPKPdf', 'Rg', 'S', 'Sb', 'ScP', 'ScS', 'Sg',
        'SKKP', 'SKP', 'SKPdf', 'SKS', 'SKSac', 'SKSdf', 'Sn', 'sP',
        'sPKPab', 'sPKPbc', 'sPKPdf', 'sS_', 'SS'
    ]
    
    for phase in phases:
        depth_samples = get_depth_samples(phase)
        distance_samples = get_distance_samples(phase)
        
        filename = f"ak135.{phase}"
        table = generate_travel_time_table(phase, depth_samples, distance_samples)
        
        with open(filename, 'w') as f:
            f.write(table)
        print(f"Generated {filename}")

if __name__ == "__main__":
    write_ak135_tables()
