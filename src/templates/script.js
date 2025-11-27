// Modern JavaScript for enhanced interactivity
document.addEventListener('DOMContentLoaded', function() {
    // Mobile navigation toggle
    const navToggle = document.getElementById('navToggle');
    const sidebar = document.getElementById('sidebar');
    
    if (navToggle && sidebar) {
        navToggle.addEventListener('click', function() {
            sidebar.classList.toggle('open');
        });
        
        // Close sidebar when clicking outside on mobile
        document.addEventListener('click', function(event) {
            if (window.innerWidth <= 768) {
                if (!sidebar.contains(event.target) && !navToggle.contains(event.target)) {
                    sidebar.classList.remove('open');
                }
            }
        });
    }
    
    // Smooth scrolling for navigation links
    const navLinks = document.querySelectorAll('.nav-link[href^="#"]');
    navLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            
            const targetId = this.getAttribute('href').substring(1);
            const targetElement = document.getElementById(targetId);
            
            if (targetElement) {
                // Update active link
                navLinks.forEach(l => l.classList.remove('active'));
                this.classList.add('active');
                
                // Smooth scroll to target
                targetElement.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
                
                // Close mobile sidebar
                if (window.innerWidth <= 768) {
                    sidebar.classList.remove('open');
                }
            }
        });
    });
    
    // Update active navigation link on scroll
    const sections = document.querySelectorAll('.content-section');
    const observerOptions = {
        root: null,
        rootMargin: '-20% 0px -70% 0px',
        threshold: 0
    };
    
    const observer = new IntersectionObserver(function(entries) {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                const targetId = entry.target.id;
                const correspondingLink = document.querySelector(`.nav-link[href="#${targetId}"]`);
                
                if (correspondingLink) {
                    navLinks.forEach(link => link.classList.remove('active'));
                    correspondingLink.classList.add('active');
                }
            }
        });
    }, observerOptions);
    
    sections.forEach(section => {
        observer.observe(section);
    });
    
    // Add loading animation for iframes
    const iframes = document.querySelectorAll('.plot-container iframe');
    iframes.forEach(iframe => {
        const container = iframe.closest('.plot-container');
        
        // Add loading indicator
        const loadingDiv = document.createElement('div');
        loadingDiv.className = 'loading-indicator';
        loadingDiv.innerHTML = `
            <div class="loading-spinner"></div>
            <p>Loading plot...</p>
        `;
        
        // Add CSS for loading indicator
        const style = document.createElement('style');
        style.textContent = `
            .loading-indicator {
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                text-align: center;
                color: var(--gray-500);
                z-index: 10;
            }
            
            .loading-spinner {
                width: 40px;
                height: 40px;
                border: 3px solid var(--gray-200);
                border-top: 3px solid var(--primary-color);
                border-radius: 50%;
                animation: spin 1s linear infinite;
                margin: 0 auto 1rem;
            }
            
            @keyframes spin {
                0% { transform: rotate(0deg); }
                100% { transform: rotate(360deg); }
            }
            
            .plot-container {
                position: relative;
            }
            
            .plot-container iframe {
                opacity: 0;
                transition: opacity 0.3s ease;
            }
            
            .plot-container iframe.loaded {
                opacity: 1;
            }
        `;
        
        if (!document.querySelector('style[data-loading-styles]')) {
            style.setAttribute('data-loading-styles', 'true');
            document.head.appendChild(style);
        }
        
        container.style.position = 'relative';
        container.appendChild(loadingDiv);
        
        iframe.addEventListener('load', function() {
            loadingDiv.remove();
            iframe.classList.add('loaded');
        });
    });
    
    // Add animation to plot cards on scroll
    const plotCards = document.querySelectorAll('.plot-card');
    const cardObserver = new IntersectionObserver(function(entries) {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.style.opacity = '1';
                entry.target.style.transform = 'translateY(0)';
            }
        });
    }, {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    });
    
    plotCards.forEach(card => {
        card.style.opacity = '0';
        card.style.transform = 'translateY(20px)';
        card.style.transition = 'opacity 0.6s ease, transform 0.6s ease';
        cardObserver.observe(card);
    });
    
    // Add keyboard navigation
    document.addEventListener('keydown', function(e) {
        // Press 'T' to toggle sidebar on mobile
        if (e.key === 't' || e.key === 'T') {
            if (window.innerWidth <= 768) {
                sidebar.classList.toggle('open');
            }
        }
        
        // Press 'Escape' to close sidebar on mobile
        if (e.key === 'Escape') {
            if (window.innerWidth <= 768) {
                sidebar.classList.remove('open');
            }
        }
    });
    
    // Add print styles
    const printStyles = document.createElement('style');
    printStyles.textContent = `
        @media print {
            .navbar,
            .sidebar,
            .footer {
                display: none !important;
            }
            
            .main-content {
                margin-left: 0 !important;
                max-width: none !important;
            }
            
            .plot-card {
                break-inside: avoid;
                page-break-inside: avoid;
            }
            
            .content-section {
                page-break-before: auto;
            }
        }
    `;
    document.head.appendChild(printStyles);
});
