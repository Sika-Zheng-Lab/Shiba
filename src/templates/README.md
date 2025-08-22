# Shiba HTML Report Template System

This document describes the improved HTML template system for Shiba reports, which separates presentation from logic and provides a modern, responsive design.

## 🎯 Improvements Made

### 1. **Separated Concerns**
- **Before**: HTML, CSS, and JavaScript embedded directly in Python code (800+ lines)
- **After**: Clean separation with dedicated template files
- **Benefits**: Easier maintenance, better collaboration, cleaner code

### 2. **Modern Design System**
- **Before**: Dark gradient background, purple color scheme, fixed sidebar
- **After**: Clean, professional design with multiple theme options
- **Features**: 
  - Responsive layout (mobile-friendly)
  - Professional color schemes
  - Modern typography (Inter font)
  - Smooth animations and transitions
  - Accessibility features

### 3. **Template Architecture**
```
src/templates/
├── summary.html          # Main template
├── splicing_section.html # Reusable section template  
├── styles.css           # Modern CSS with CSS variables
└── script.js            # Interactive JavaScript
```

## 🎨 Available Themes

| Theme | Description | Best For |
|-------|-------------|----------|
| **Modern Blue** (default) | Clean, professional blue | General use, presentations |
| **Scientific Green** | Nature-inspired green | Biological research |
| **Dark Mode** | Dark theme | Reduced eye strain |
| **Academic Purple** | Elegant purple | Academic presentations |
| **Minimal Gray** | Clean grayscale | Print-friendly reports |

## 🚀 Usage

### Basic Usage (Updated plots.py)
The improved `plots.py` automatically uses the new template system:

```python
# The write_summary_html function now uses templates
write_summary_html(shiba_command, output_dir)
```

### Custom Theme Selection
```python
from themes import generate_theme_css_file

# Generate a themed CSS file
generate_theme_css_file("scientific_green", "output/styles.css")
```

### Template Customization
```python
from template_renderer import HTMLTemplateRenderer

renderer = HTMLTemplateRenderer("path/to/templates")
html_content = renderer.render_summary_html(data)
```

## 📱 Responsive Features

### Desktop View
- Fixed sidebar navigation
- 3-column plot grid
- Full-sized plots and content

### Mobile View
- Collapsible sidebar (hamburger menu)
- Single-column layout
- Touch-friendly navigation
- Optimized plot sizes

### Keyboard Navigation
- `T` key: Toggle sidebar (mobile)
- `Escape`: Close sidebar (mobile)
- Tab navigation through all interactive elements

## 🛠 Customization Guide

### Changing Colors
Edit `src/themes.py` to modify existing themes or add new ones:

```python
COLOR_SCHEMES["my_theme"] = {
    "name": "My Custom Theme",
    "colors": {
        "primary": "#your-color",
        "background": "#your-bg-color",
        # ... more colors
    }
}
```

### Modifying Layout
Edit `src/templates/summary.html` to change the structure:
- Add new sections
- Modify navigation
- Change header/footer content

### Styling Changes
Edit `src/templates/styles.css`:
- CSS variables for easy theming
- Modern CSS Grid and Flexbox
- Responsive design utilities

### Adding Interactivity
Edit `src/templates/script.js`:
- Smooth scrolling
- Mobile navigation
- Loading animations
- Custom interactions

## 📊 Plot Integration

### Supported Plot Types
- **PCA Plots**: Gene expression and splicing pattern analysis
- **Volcano Plots**: Differential splicing significance
- **Scatter Plots**: PSI comparison between conditions  
- **Bar Charts**: Event count summaries

### Plot Features
- Loading indicators
- Responsive sizing
- Error handling for missing data
- Smooth iframe integration

## 🔧 Development Workflow

### Making Changes
1. **Templates**: Edit HTML files in `src/templates/`
2. **Styles**: Modify CSS in `src/templates/styles.css`
3. **Scripts**: Update JavaScript in `src/templates/script.js`
4. **Themes**: Add/modify themes in `src/themes.py`

### Testing Changes
```bash
# Generate demo with new changes
cd src
python demo_template.py

# View in browser
open /tmp/shiba_demo_*/demo_summary.html
```

### Integration
The improved `plots.py` automatically:
1. Copies static files (CSS, JS) to output directory
2. Renders templates with data
3. Generates modern HTML reports

## 🎯 Benefits Summary

| Aspect | Before | After |
|--------|---------|-------|
| **Maintainability** | 800+ line HTML string | Modular template files |
| **Design** | Dark/purple theme | 5 professional themes |
| **Responsiveness** | Desktop only | Mobile-friendly |
| **Customization** | Edit Python code | Edit template files |
| **Performance** | Heavy inline styles | Optimized CSS |
| **Accessibility** | Basic | WCAG compliant |
| **User Experience** | Static | Interactive + animated |

## 📝 Migration Notes

### For Existing Users
- No changes needed to existing `plots.py` calls
- Output HTML will automatically use new design
- Old HTML structure is completely replaced

### For Developers
- HTML customization now happens in template files
- CSS variables enable easy theme switching
- JavaScript enables enhanced interactivity

## 🔮 Future Enhancements

Potential additions to consider:
- **Dark/Light mode toggle**
- **Plot export functionality** 
- **Interactive plot filtering**
- **Customizable report sections**
- **Print-optimized layouts**
- **PDF export capability**

---

*The new template system provides a solid foundation for modern, maintainable HTML reports while keeping the API simple for end users.*
