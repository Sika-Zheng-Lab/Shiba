"""
Color schemes and theme configurations for Shiba HTML reports.
You can easily switch between different themes by changing the ACTIVE_THEME variable.
"""

# Available color schemes
COLOR_SCHEMES = {
    "modern_blue": {
        "name": "Modern Blue",
        "description": "Clean, professional blue theme",
        "colors": {
            "primary": "#2563eb",
            "primary_dark": "#1d4ed8", 
            "secondary": "#64748b",
            "accent": "#0ea5e9",
            "success": "#10b981",
            "warning": "#f59e0b",
            "error": "#ef4444",
            "background": "#f8fafc",
            "surface": "#ffffff",
            "text": "#334155",
            "text_light": "#64748b"
        }
    },
    
    "scientific_green": {
        "name": "Scientific Green", 
        "description": "Nature-inspired green theme for biological data",
        "colors": {
            "primary": "#059669",
            "primary_dark": "#047857",
            "secondary": "#6b7280", 
            "accent": "#10b981",
            "success": "#22c55e",
            "warning": "#f59e0b",
            "error": "#ef4444",
            "background": "#f0fdf4",
            "surface": "#ffffff",
            "text": "#1f2937",
            "text_light": "#6b7280"
        }
    },
    
    "dark_mode": {
        "name": "Dark Mode",
        "description": "Dark theme for reduced eye strain",
        "colors": {
            "primary": "#3b82f6",
            "primary_dark": "#2563eb",
            "secondary": "#9ca3af",
            "accent": "#06b6d4", 
            "success": "#10b981",
            "warning": "#f59e0b",
            "error": "#ef4444",
            "background": "#111827",
            "surface": "#1f2937",
            "text": "#f9fafb",
            "text_light": "#d1d5db"
        }
    },
    
    "academic_purple": {
        "name": "Academic Purple",
        "description": "Elegant purple theme for academic presentations",
        "colors": {
            "primary": "#7c3aed",
            "primary_dark": "#6d28d9",
            "secondary": "#64748b",
            "accent": "#8b5cf6",
            "success": "#10b981", 
            "warning": "#f59e0b",
            "error": "#ef4444",
            "background": "#faf5ff",
            "surface": "#ffffff",
            "text": "#1e1b4b",
            "text_light": "#64748b"
        }
    },
    
    "minimal_gray": {
        "name": "Minimal Gray",
        "description": "Clean, minimal grayscale theme",
        "colors": {
            "primary": "#374151",
            "primary_dark": "#1f2937",
            "secondary": "#9ca3af",
            "accent": "#6b7280",
            "success": "#10b981",
            "warning": "#f59e0b", 
            "error": "#ef4444",
            "background": "#f9fafb",
            "surface": "#ffffff",
            "text": "#111827",
            "text_light": "#6b7280"
        }
    }
}

# Set the active theme here
ACTIVE_THEME = "modern_blue"

def get_active_theme():
    """Get the currently active color scheme."""
    return COLOR_SCHEMES[ACTIVE_THEME]

def get_css_variables(theme_name=None):
    """Generate CSS custom properties for the specified theme."""
    theme = COLOR_SCHEMES.get(theme_name, COLOR_SCHEMES[ACTIVE_THEME])
    colors = theme["colors"]
    
    css_vars = ":root {\n"
    css_vars += f"    /* {theme['name']} - {theme['description']} */\n"
    
    for color_name, color_value in colors.items():
        css_var_name = color_name.replace("_", "-")
        css_vars += f"    --{css_var_name}: {color_value};\n"
    
    css_vars += "}\n"
    return css_vars

def generate_theme_css_file(theme_name, output_path):
    """Generate a complete CSS file for a specific theme."""
    
    # Read the base CSS template
    base_css_path = "/rhome/naotok/Shiba/src/templates/styles.css"
    with open(base_css_path, 'r') as f:
        base_css = f.read()
    
    # Get theme variables
    theme_vars = get_css_variables(theme_name)
    
    # Replace the :root section in base CSS
    import re
    pattern = r':root\s*{[^}]*}'
    updated_css = re.sub(pattern, theme_vars.strip(), base_css, count=1)
    
    # Write the themed CSS
    with open(output_path, 'w') as f:
        f.write(updated_css)
    
    return output_path

def list_available_themes():
    """List all available themes with descriptions."""
    print("Available Shiba HTML Report Themes:")
    print("=" * 40)
    
    for theme_key, theme_data in COLOR_SCHEMES.items():
        status = " (ACTIVE)" if theme_key == ACTIVE_THEME else ""
        print(f"• {theme_data['name']}{status}")
        print(f"  Key: {theme_key}")
        print(f"  Description: {theme_data['description']}")
        print()

if __name__ == "__main__":
    list_available_themes()
    
    # Example of generating theme files
    import os
    theme_dir = "/tmp/shiba_themes"
    os.makedirs(theme_dir, exist_ok=True)
    
    print(f"Generating theme CSS files in {theme_dir}:")
    for theme_key in COLOR_SCHEMES.keys():
        output_file = os.path.join(theme_dir, f"styles_{theme_key}.css")
        generate_theme_css_file(theme_key, output_file)
        print(f"  ✓ {theme_key} -> {output_file}")
